"""Stand-alone utilities for computing and plotting Keplerian orbits."""
import numpy as np

try:
    from numba import njit

    do_numba = True
except ModuleNotFoundError:
    njit = None
    do_numba = False

pi = np.pi
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)


def decorator_if(decorator, condition):
    """Decorate a function with a decorator if a condition is satisfied.

    Otherwise leave the function undecorated.
    """

    def dec(func):
        if condition:
            return decorator(func)
        else:
            return func

    return dec


def njit_if_numba():
    """Decorate a function with njit if available."""
    return decorator_if(njit, do_numba)


@njit_if_numba()
def root_sum_sq(A1_x):
    """Return the root of the summed squares of all x in A1_x (incl. 1D x)."""
    ans = A1_x[0] ** 2
    for x in A1_x[1:]:
        ans += x**2

    return np.sqrt(ans)


@njit_if_numba()
def mod_angle(a):
    """Set an angle between 0 and 2pi"""
    return np.mod(2 * pi + np.mod(a, 2 * pi), 2 * pi)


@njit_if_numba()
def arccos2(o, a, disambiguator=1):
    """Return arccos(o/a), using the disambiguator to set the quadrant."""
    cosine = o / a
    # In normal bounds
    if cosine > -1 and cosine < 1:
        angle = np.arccos(cosine)
        if disambiguator < 0:
            angle = -angle
    # Exceptions
    elif cosine <= -1:
        angle = pi
    else:
        angle = 0

    return angle


@njit_if_numba()
def _orbit_from_pos_vel_M_p_m(A1_pos, A1_vel, M_p, m, tol=1e-8):
    """Calculate orbital parameters from the position and velocity.

    Parameters
    ----------
    A1_pos : [float]
        The x, y, z coordinates relative to the primary (m).

    A1_vel : [float]
        The x, y, z velocities relative to the primary (m s^-1).

    M_p : float
        The mass of the primary (kg).

    m : float
        The mass of the orbiting body (kg).

    tol : float (opt.)
        Tolerance for parabolic e and planar i values.

    Returns
    -------
    r, v, v_r, a, h, e, n, i, nu, E, M, theta, Omega, pomega, omega : float
        See definitions in the Orbit class.
    """
    mu = G * (M_p + m)

    # Distance and speed
    r = root_sum_sq(A1_pos)
    v = root_sum_sq(A1_vel)

    # Radial velocity
    v_r = np.vdot(A1_pos, A1_vel) / r

    # Specific angular momentum
    A1_h = np.cross(A1_pos, A1_vel)
    h = root_sum_sq(A1_h)

    # Eccentricity (vector from apoapsis to periapsis)
    A1_e = np.cross(A1_vel, A1_h) / mu - A1_pos / r
    e = root_sum_sq(A1_e)

    # Inclination
    i = np.arccos(A1_h[2] / h)

    # Angles
    # Treat ~zero inclination separately for the well or not-well defined angles
    if i < tol or i > pi - tol:
        # Roughly zero inclination, so use longitudes instead of node angles
        Omega = 0
        # True longitude
        theta = arccos2(A1_pos[0], r, A1_pos[1])
        # Longitude of periapsis
        pomega = arccos2(A1_e[0], e, A1_e[1])
        # Argument of periapsis
        omega = pomega - Omega
        # True anomaly
        nu = theta - pomega

        # Retrograde changes signs
        if i > pi / 2:
            omega = -omega
            nu = -nu
    else:
        # Vector pointing to the ascending node
        A1_n = np.array([-A1_h[1], A1_h[0], 0])
        n_ = root_sum_sq(A1_n)
        # Longitude of ascending node
        # If y component < 0, then Omega > pi
        Omega = arccos2(-A1_h[1], n_, A1_n[1])

        # Inclined, use the combined angles well-defined in the orbital plane
        omega_plus_nu = arccos2(np.vdot(A1_n, A1_pos), n_ * r, A1_pos[2])
        # Argument of periapsis
        omega = arccos2(np.vdot(A1_n, A1_e), n_ * e, A1_e[2])
        # True anomaly
        nu = omega_plus_nu - omega

        # Retrograde changes signs
        if i < pi / 2:
            sign = 1
        else:
            sign = -1
        # Longitude of periapsis
        pomega = Omega + sign * omega
        # True longitude
        theta = Omega + sign * omega_plus_nu

    # Eccentric anomaly and mean anomaly
    # Treat parabolic separately since semi-major axis not defined
    if abs(e - 1) < tol:
        a = np.inf

        D = np.tan(nu / 2)
        if v_r < 0:
            D = 2 * pi - D
        M = D + D**3 / 3

        # Label D as eccentic anomaly for returned values
        E = D
    else:
        # Semi-major axis
        a = 1 / (2 / r - v**2 / mu)

        # Mean motion (negative if hyperbolic)
        n = np.sign(a) * np.sqrt(mu / abs(a) ** 3)

        # Eccentric anomaly and mean anomaly
        if e < 1:
            # If v_r < 0, then after apoapsis, so E > pi
            E = arccos2((1 - r / a) / e, v_r)
            M = E - e * np.sin(E)
        else:
            E = np.arccosh((1 - r / a) / e)
            if v_r < 0:
                E = 2 * pi - E
            M = e * np.sinh(E) - E

    return (
        r,
        v,
        v_r,
        a,
        h,
        e,
        n,
        i,
        mod_angle(nu),
        mod_angle(E),
        mod_angle(M),
        mod_angle(theta),
        mod_angle(Omega),
        mod_angle(pomega),
        mod_angle(omega),
    )


class Orbit:
    """Orbital parameters for a body orbiting a primary.

    Provide A1_pos, A1_vel, M_p, m to automatically derive the other parameters.

    Parameters (all optional)
    ----------
    A1_pos : [float]
        The x, y, z coordinates relative to the primary (m).

    A1_vel : [float]
        The x, y, z velocities relative to the primary (m s^-1).

    M_p : float
        The mass of the primary (kg).

    m : float
        The mass of the orbiting body (kg).

    a : float
        The semi-major axis, negative for hyperbolic orbits (m).

    h : float
        The specific angular momentum (m^2 s^-1).

    e : float
        The eccentricity.

    i : float
        The inclination (rad): angle of the angular momentum from the z axis.

    nu : float
        The true anomaly (rad): angle away from periapsis.

    theta : float
        The true longitude (rad): angle from x axis if inclination were zero.

    Omega : float
        The longitude of the ascending node (rad): angle from the x axis of
        where the orbit passes up through the z=0 plane.

    pomega : float
        The longitude of periapsis (rad): angle of periapsis from the x axis.

    omega : float
        The argument of periapsis (rad): angle of periapsis from the ascending
        node.

    q : float
        The periapsis (m), which for parabolic orbits is the defining parameter
        instead of a.

    id : int
        An ID number.

    Other attributes (see also properties)
    ----------------
    r : float
        Distance from primary (m).

    v : float
        The speed (m s^-1).

    v_r : float
        The radial speed (m s^-1).

    n : float
        The mean motion (rad s^-1): average angular speed.

    M : float
        The mean anomaly (rad): fraction of period elapsed since periapsis.

    E : float
        The eccentric anomaly (rad): angle from the ellipse centre.
    """

    def __init__(
        self,
        A1_pos=None,
        A1_vel=None,
        M_p=None,
        m=None,
        a=None,
        h=None,
        e=None,
        i=None,
        nu=None,
        theta=None,
        Omega=None,
        pomega=None,
        omega=None,
        q=None,
        id=None,
    ):
        self.A1_pos = A1_pos
        self.A1_vel = A1_vel
        self.M_p = M_p
        self.m = m
        self.a = a
        self.h = h
        self.e = e
        self.i = i
        self.nu = nu
        self.theta = theta
        self.Omega = Omega
        self.pomega = pomega
        self.omega = omega
        self._q = q
        self.id = id

        # Tolerance for planar, circular, and parabolic orbits
        self.tol = 1e-8

        # Derive other parameters if the relevant inputs were provided
        if (
            A1_pos is not None
            and A1_vel is not None
            and M_p is not None
            and m is not None
            and a is None
            and e is None
        ):
            self.orbit_from_pos_vel_M_p_m(A1_pos, A1_vel, M_p, m)

    def orbit_from_pos_vel_M_p_m(self, A1_pos, A1_vel, M_p, m):
        """See _orbit_from_pos_vel_M_p_m()."""

        if any(np.isnan(A1_pos)) or any(np.isnan(A1_vel)):
            return self.nan_orbit()

        # Ensure types for numba
        self.A1_pos = np.array(A1_pos, dtype=np.float64)
        self.A1_vel = np.array(A1_vel, dtype=np.float64)
        self.M_p = M_p
        self.m = m

        # Zeros break these equations
        assert not all(self.A1_pos == 0)
        assert not all(self.A1_vel == 0)

        # Compute and set attributes
        (
            self.r,
            self.v,
            self.v_r,
            self.a,
            self.h,
            self.e,
            self.n,
            self.i,
            self.nu,
            self.E,
            self.M,
            self.theta,
            self.Omega,
            self.pomega,
            self.omega,
        ) = _orbit_from_pos_vel_M_p_m(
            self.A1_pos, self.A1_vel, self.M_p, self.m, tol=self.tol
        )

        return self

    def nan_orbit(self):
        """Set all attributes to nan"""
        self.A1_pos = np.nan
        self.A1_vel = np.nan
        self.M_p = np.nan
        self.m = np.nan
        self.r = np.nan
        self.v = np.nan
        self.v_r = np.nan
        self.a = np.nan
        self.h = np.nan
        self.e = np.nan
        self.n = np.nan
        self.i = np.nan
        self.nu = np.nan
        self.E = np.nan
        self.M = np.nan
        self.theta = np.nan
        self.Omega = np.nan
        self.pomega = np.nan
        self.omega = np.nan

        return self

    def t_from_nu(self, nu):
        """Calculate the time since periapsis from the true anomaly."""
        # Mean anomaly
        M = self.M_from_nu(nu)

        # Time
        if self.is_parabolic:
            t = M * np.sqrt(2 * self.q**3 / self.mu)
        else:
            t = M * np.sqrt(abs(self.a**3) / self.mu)

        return t

    def nu_from_t(self, t):
        """Calculate the true anomaly from the time since periapsis."""
        # Mean anomaly
        if self.is_parabolic:
            M = t * np.sqrt(self.mu / (2 * self.q**3))
        else:
            M = t * np.sqrt(self.mu / abs(self.a**3))

        # True anomaly
        nu = self.nu_from_M(M)

        return mod_angle(nu)

    def pos_vel_from_nu(self, nu, do_vel=True):
        """Calculate position and velocity from the orbit and true anomaly.

        Parameters
        ----------
        nu : float
            The true anomaly (rad): angle away from periapsis.

        do_vel : bool (opt.)
            If False, then don't return the velocity. Removes the requirement of
            M_p, m attributes.

        Required attributes
        -------------------
        M_p, m, a, e, i, Omega, omega (see class docstring), or q instead
        of a for parabolic orbits.

        Returns
        -------
        A1_pos, A1_vel : [float]
            The position (m) and velocity (m s^-1).
        """
        # Modified from rebound reb_tools_orbit_to_particle_err().
        # Distance
        if abs(self.e - 1) < self.tol:
            # Parabolic
            r = 2 * self.q / (1 + np.cos(nu))

            if do_vel:
                v0 = np.sqrt(2 * self.mu / r)
        else:
            r = self.a * (1 - self.e**2) / (1 + self.e * np.cos(nu))

            # This form works for elliptical and hyperbolic orbits
            if do_vel:
                v0 = np.sqrt(self.mu / (self.a * (1 - self.e**2)))

        cn = np.cos(nu)
        sn = np.sin(nu)
        ci = np.cos(self.i)
        si = np.sin(self.i)
        cO = np.cos(self.Omega)
        sO = np.sin(self.Omega)
        co = np.cos(self.omega)
        so = np.sin(self.omega)

        A1_pos = np.zeros(3)
        A1_vel = np.zeros(3)

        # Murray & Dermott Eq 2.122
        A1_pos[0] = r * (cO * (co * cn - so * sn) - sO * (so * cn + co * sn) * ci)
        A1_pos[1] = r * (sO * (co * cn - so * sn) + cO * (so * cn + co * sn) * ci)
        A1_pos[2] = r * (so * cn + co * sn) * si

        if do_vel:
            # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from
            # S2.8 to the velocities in the orbital plane
            A1_vel[0] = v0 * (
                (self.e + cn) * (-ci * co * sO - cO * so)
                - sn * (co * cO - ci * so * sO)
            )
            A1_vel[1] = v0 * (
                (self.e + cn) * (ci * co * cO - sO * so) - sn * (co * sO + ci * so * cO)
            )
            A1_vel[2] = v0 * ((self.e + cn) * co * si - sn * si * so)

            return A1_pos, A1_vel
        else:
            return A1_pos

    def pos_vel_from_t(self, t, do_vel=True):
        """Calculate position and velocity from the orbit and time since periapsis."""
        return self.pos_vel_from_nu(nu=self.nu_from_t(t), do_vel=do_vel)

    def E_from_nu(self, nu):
        """Calculate the eccentric anomaly from the true anomaly."""
        if self.e < 1:
            E = np.arccos((self.e + np.cos(nu)) / (1 + self.e * np.cos(nu)))
        else:
            E = np.arccosh((self.e + np.cos(nu)) / (1 + self.e * np.cos(nu)))

        return mod_angle(E)

    def E_from_M(self, M):
        """Calculate the eccentric anomaly from the mean anomaly."""
        M = mod_angle(M)
        e = self.e
        if e < 1:
            # Starting value
            if e < 0.8:
                E = M
            else:
                E = pi
            F = E - e * np.sin(E) - M

            # Iterate
            for iter in range(100):
                E = E - F / (1 - e * np.cos(E))
                F = E - e * np.sin(E) - M
                if abs(F) < 1e-16:
                    break
        else:
            # Starting value
            E = np.sign(M) * np.log(2 * abs(M) / e + 1.8)
            F = E - e * np.sinh(E) + M

            # Iterate
            for iter in range(100):
                E = E - F / (1 - e * np.cosh(E))
                F = E - e * np.sinh(E) + M
                if abs(F) < 1e-16:
                    break

        return mod_angle(E)

    def M_from_E(self, E):
        """Calculate the mean anomaly from the eccentric anomaly."""
        if self.e < 1:
            return E - self.e * np.sin(E)
        else:
            return self.e * np.sinh(E) - E

        return mod_angle(nu)

    def M_from_nu(self, nu):
        """Calculate the mean anomaly from the true anomaly."""
        E = self.E_from_nu(nu)

        return self.M_from_E(E)

    def nu_from_E(self, E):
        """Calculate the true anomaly from the eccentric anomaly."""
        if self.e < 1:
            nu = 2 * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E / 2))
        else:
            nu = 2 * np.arctan(np.sqrt((1 + self.e) / (self.e - 1)) * np.tanh(E / 2))

        return mod_angle(nu)

    def nu_from_M(self, M):
        """Calculate the true anomaly from the mean anomaly,."""
        E = self.E_from_M(M)

        return self.nu_from_E(E)

    def nu_from_r(self, r):
        """Calculate the true anomaly from the radial distance."""
        if self.q > r or self.Q < r:
            return np.nan
        elif self.is_parabolic:
            return np.arccos(2 * self.q / r - 1)
        else:
            return np.arccos((self.a * (1 - self.e**2) / r - 1) / self.e)

    def r_from_nu(self, nu):
        """Calculate the radial distance from the true anomaly."""
        if self.is_parabolic:
            return 2 * self.q / (1 + np.cos(nu))
        else:
            return self.a * (1 - self.e**2) / (1 + self.e * np.cos(nu))

    def t_from_r(self, r):
        """Calculate the time since periapsis from the radial distance."""
        return self.t_from_nu(self.nu_from_r(r))

    def r_from_t(self, t):
        """Calculate the radial distance from the time since periapsis."""
        return self.r_from_nu(self.nu_from_t(t))

    def param(self, param):
        """Return a parameter selected by its name."""
        return getattr(self, param)

    @property
    def mu(self):
        """Gravitational parameter (m^3 s^-2)."""
        return G * (self.M_p + self.m)

    @property
    def A1_pos_bary(self):
        """Position of the barycentre relative to the primary (m)."""
        return self.A1_pos * self.m / (self.M_p + self.m)

    @property
    def A1_pos_q(self):
        """Position of the periapsis relative to the primary (m)."""
        return self.pos_vel_from_nu(nu=0, do_vel=False)

    @property
    def A1_pos_Q(self):
        """Position of the apoapsis relative to the primary (m)."""
        return self.pos_vel_from_nu(nu=pi, do_vel=False)

    @property
    def period(self):
        """Period (s)."""
        if self.e < 1:
            return np.sqrt(4 * pi**2 * self.a**3 / self.mu)
        else:
            return np.nan

    @property
    def is_parabolic(self):
        """Is it (close enough to) parabolic?"""
        try:
            v_inf_2 = self.v**2 - 2 * self.mu / self.r
            return np.round(abs(v_inf_2), 7) == 0
        except AttributeError:
            return abs(self.e - 1) < self.tol

    @property
    def b(self):
        """Semi-minor axis (m)."""
        if self.is_parabolic:
            return np.nan
        else:
            if self.e < 1:
                return self.a * np.sqrt(1 - self.e**2)
            else:
                return self.a * np.sqrt(self.e**2 - 1)

    @property
    def p(self):
        """Semi-latus rectum (m)."""
        return self.h**2 / self.mu

    @property
    def q(self):
        """Periapsis (m)."""
        if hasattr(self, "_q") and self._q is not None and not np.isnan(self._q):
            return self._q
        if self.is_parabolic:
            return self.p / 2
        else:
            return self.a * (1 - self.e)

    @property
    def Q(self):
        """Apoapsis (m)."""
        if self.is_parabolic:
            return np.inf
        else:
            return self.a * (1 + self.e)

    @property
    def R_Hill(self):
        """Hill radius (m)."""
        return self.a * np.cbrt(self.m / (3 * self.M_p))

    @property
    def v_esc_p(self):
        """Escape speed at periapsis (m s^-1)."""
        return np.sqrt(2 * self.mu / self.q)

    @property
    def v_p(self):
        """Speed at periapsis (m s^-1)."""
        if self.is_parabolic:
            return self.v_esc_p
        else:
            return np.sqrt(self.mu * (2 / self.q - 1 / self.a))

    @property
    def eps(self):
        """Specific energy (J kg^-1)."""
        if self.is_parabolic:
            return 0
        else:
            return -self.mu / (2 * self.a)

    @property
    def v_inf(self):
        """Speed at infinity (m s^-1)."""
        if self.e > 1:
            return np.sqrt(self.v**2 - 2 * self.mu / self.r)
        else:
            return np.nan

    @property
    def l(self):
        """Mean longitude (rad)."""
        # Treat zero inclination separately given the well or not-well defined angles
        if self.i < self.tol or self.i > pi - self.tol:
            # Prograde or retrograde changes signs
            if self.i < pi / 2:
                # Prograde
                if self.e > self.tol:
                    l = self.pomega + self.M
                else:
                    # pomega not well defined for circular orbits
                    l = self.theta - 2 * self.e * np.sin(self.nu)
            else:
                # Retrograde
                if self.e > self.tol:
                    l = self.pomega - self.M
                else:
                    # pomega not well defined for circular orbits
                    l = self.theta + 2 * self.e * np.sin(self.nu)
        else:
            # Prograde or retrograde changes signs
            if self.i < pi / 2:
                # Prograde
                if self.e > self.tol:
                    l = self.pomega + self.M
                else:
                    # pomega not well defined for circular orbits
                    l = self.theta - 2 * self.e * np.sin(self.nu)
            else:
                # Retrograde
                if self.e > self.tol:
                    l = self.pomega - self.M
                else:
                    # pomega not well defined for circular orbits
                    l = self.theta + 2 * self.e * np.sin(self.nu)

        return mod_angle(l)

    @property
    def alpha(self):
        """Angle of velocity away from radial (rad)."""
        return np.arccos(np.dot(self.A1_vel, self.A1_pos) / (self.v * self.r))

    @property
    def r_mean(self):
        """Mean distance (m), averaged over time."""
        if self.e < 1:
            return self.a * (1 + self.e**2 / 2)
        else:
            return np.inf

    @property
    def perimeter(self):
        """Perimeter of the orbital ellipse (m)."""
        if self.e < 1:
            # Infinite series approximation
            return (
                2
                * self.a
                * pi
                * (
                    1
                    - 1 / 4 * self.e**2
                    - 3 / 64 * self.e**4
                    - 5 / 256 * self.e**6
                )
            )
        else:
            return np.inf

    @property
    def v_mean(self):
        """Mean speed (m s^-1), averaged over time."""
        if self.e < 1:
            return self.perimeter / self.period
        else:
            return self.v_inf


class OrbitSet:
    """A set of Orbit objects.

    Parameters
    ----------
    A1_o or A2_o : [Orbit] or [[Orbit]] (opt.)
        The 1D or 2D array of Orbit objects. Both are set as attributes, either
        by adding a dimension to A1_o or by flattening A2_o.

    <or>
    A2_pos, A2_vel, A1_M_p, A1_m : [[float]], [[float]], [float], [float] (opt.)
    A3_pos, A3_vel, A2_M_p, A2_m : [[[float]]], [[[float]]], [[float]], [[float]] (opt.)
        1D or 2D arrays of the four parameters for orbit_from_pos_vel_M_p_m(),
        used to create the corresponding A1_o or A2_o array of Orbit objects.

    A1_id or A2_id : [int] or [[int]] (opt.)
        1D or 2D arrays of ID numbers for the corresponding Orbit objects.

    Attributes
    ----------
    num_orb : int
        The total number of Orbit objects.

    num_orb_1, num_orb_2 : int
        The dimensions of the 2D Orbit array.

    Methods
    -------
    A1_*, A2_* : [float], [[float]]
        Return a convenient array of parameters for each Orbit, with the same
        dimensions as the corresponding A1_o or A2_o.
    """

    def __init__(
        self,
        A1_o=None,
        A2_o=None,
        A2_pos=None,
        A2_vel=None,
        A1_M_p=None,
        A1_m=None,
        A1_id=None,
        A3_pos=None,
        A3_vel=None,
        A2_M_p=None,
        A2_m=None,
        A2_id=None,
    ):
        # Create the Orbit objects if not provided
        if A1_o is None and A2_o is None:
            # 1D input
            if A1_m is not None:
                A1_o = [
                    Orbit(
                        A1_pos=A2_pos[i],
                        A1_vel=A2_vel[i],
                        M_p=A1_M_p[i],
                        m=A1_m[i],
                    )
                    for i in range(len(A1_m))
                ]
            # 2D input
            else:
                assert A2_m is not None
                A2_o = [
                    [
                        Orbit(
                            A1_pos=A3_pos[i, j],
                            A1_vel=A3_vel[i, j],
                            M_p=A2_M_p[i, j],
                            m=A2_m[i, j],
                        )
                        for j in range(len(A2_m[0]))
                    ]
                    for i in range(len(A2_m))
                ]

        # Parse the Orbit objects
        # 1D input
        if A1_o is not None:
            A1_o = np.array(A1_o, dtype=object)

            # Set IDs if provided
            if A1_id is not None:
                for i in range(len(A1_o)):
                    A1_o[i].id = A1_id[i]

            self.A1_o = A1_o
            self.A2_o = np.array([self.A1_o], dtype=object)
        # 2D input
        else:
            assert A2_o is not None
            A2_o = np.array(A2_o, dtype=object)

            # Set IDs if provided
            if A1_id is not None:
                for i in range(len(A2_o)):
                    for j in range(len(A2_o[0])):
                        A2_o[i, j].id = A1_id[i, j]

            self.A2_o = A2_o
            self.A1_o = np.array([o for A1_o in self.A2_o for o in A1_o], dtype=object)

        # Sanity checks
        assert isinstance(self.A1_o[0], Orbit)
        assert isinstance(self.A2_o[0, 0], Orbit)
        self.num_orb = len(self.A1_o)
        self.num_orb_1 = len(self.A2_o)
        self.num_orb_2 = len(self.A2_o[0])
        assert self.num_orb == self.num_orb_1 * self.num_orb_2

    def A1_param(self, param):
        return np.array([o.param(param) for o in self.A1_o])

    def A2_param(self, param):
        return np.array([[o.param(param) for o in A1_o] for A1_o in self.A2_o])

    @property
    def A2_pos(self):
        return np.array([o.A1_pos for o in self.A1_o])

    @property
    def A3_pos(self):
        return np.array([[o.A1_pos for o in A1_o] for A1_o in self.A2_o])

    @property
    def A2_vel(self):
        return np.array([o.A1_vel for o in self.A1_o])

    @property
    def A3_vel(self):
        return np.array([[o.A1_vel for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_M_p(self):
        return np.array([o.M_p for o in self.A1_o])

    @property
    def A2_M_p(self):
        return np.array([[o.M_p for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_m(self):
        return np.array([o.m for o in self.A1_o])

    @property
    def A2_m(self):
        return np.array([[o.m for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_a(self):
        return np.array([o.a for o in self.A1_o])

    @property
    def A2_a(self):
        return np.array([[o.a for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_h(self):
        return np.array([o.h for o in self.A1_o])

    @property
    def A2_h(self):
        return np.array([[o.h for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_e(self):
        return np.array([o.e for o in self.A1_o])

    @property
    def A2_e(self):
        return np.array([[o.e for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_i(self):
        return np.array([o.i for o in self.A1_o])

    @property
    def A2_i(self):
        return np.array([[o.i for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_nu(self):
        return np.array([o.nu for o in self.A1_o])

    @property
    def A2_nu(self):
        return np.array([[o.nu for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_theta(self):
        return np.array([o.theta for o in self.A1_o])

    @property
    def A2_theta(self):
        return np.array([[o.theta for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_Omega(self):
        return np.array([o.Omega for o in self.A1_o])

    @property
    def A2_Omega(self):
        return np.array([[o.Omega for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_pomega(self):
        return np.array([o.pomega for o in self.A1_o])

    @property
    def A2_pomega(self):
        return np.array([[o.pomega for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_omega(self):
        return np.array([o.omega for o in self.A1_o])

    @property
    def A2_omega(self):
        return np.array([[o.omega for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_r(self):
        return np.array([o.r for o in self.A1_o])

    @property
    def A2_r(self):
        return np.array([[o.r for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_v(self):
        return np.array([o.v for o in self.A1_o])

    @property
    def A2_v(self):
        return np.array([[o.v for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_v_r(self):
        return np.array([o.v_r for o in self.A1_o])

    @property
    def A2_v_r(self):
        return np.array([[o.v_r for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_n(self):
        return np.array([o.n for o in self.A1_o])

    @property
    def A2_n(self):
        return np.array([[o.n for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_M(self):
        return np.array([o.M for o in self.A1_o])

    @property
    def A2_M(self):
        return np.array([[o.M for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_E(self):
        return np.array([o.E for o in self.A1_o])

    @property
    def A2_E(self):
        return np.array([[o.E for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_mu(self):
        return np.array([o.mu for o in self.A1_o])

    @property
    def A2_mu(self):
        return np.array([[o.mu for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_period(self):
        return np.array([o.period for o in self.A1_o])

    @property
    def A2_period(self):
        return np.array([[o.period for o in A1_o] for A1_o in self.A2_o])

    @property
    def A2_pos_q(self):
        return np.array([o.A1_pos_q for o in self.A1_o])

    @property
    def A3_pos_q(self):
        return np.array([[o.A1_pos_q for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_b(self):
        return np.array([o.b for o in self.A1_o])

    @property
    def A2_b(self):
        return np.array([[o.b for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_p(self):
        return np.array([o.p for o in self.A1_o])

    @property
    def A2_p(self):
        return np.array([[o.p for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_q(self):
        return np.array([o.q for o in self.A1_o])

    @property
    def A2_q(self):
        return np.array([[o.q for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_Q(self):
        return np.array([o.Q for o in self.A1_o])

    @property
    def A2_Q(self):
        return np.array([[o.Q for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_R_Hill(self):
        return np.array([o.R_Hill for o in self.A1_o])

    @property
    def A2_R_Hill(self):
        return np.array([[o.R_Hill for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_v_esc_p(self):
        return np.array([o.v_esc_p for o in self.A1_o])

    @property
    def A2_v_esc_p(self):
        return np.array([[o.v_esc_p for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_v_p(self):
        return np.array([o.v_p for o in self.A1_o])

    @property
    def A2_v_p(self):
        return np.array([[o.v_p for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_eps(self):
        return np.array([o.eps for o in self.A1_o])

    @property
    def A2_eps(self):
        return np.array([[o.eps for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_v_inf(self):
        return np.array([o.v_inf for o in self.A1_o])

    @property
    def A2_v_inf(self):
        return np.array([[o.v_inf for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_l(self):
        return np.array([o.l for o in self.A1_o])

    @property
    def A2_l(self):
        return np.array([[o.l for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_alpha(self):
        return np.array([o.alpha for o in self.A1_o])

    @property
    def A2_alpha(self):
        return np.array([[o.alpha for o in A1_o] for A1_o in self.A2_o])

    @property
    def A1_id(self):
        return np.array([o.id for o in self.A1_o])

    @property
    def A2_id(self):
        return np.array([[o.id for o in A1_o] for A1_o in self.A2_o])


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  orbits.py  ====\n")
