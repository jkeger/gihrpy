"""
General utilities for GIHR projects.
"""

from jkeger import *


# ========
# Initialise parameters for global access (e.g. gihr.py flag arguments)
# ========
do_silent = False
do_paper = False
do_skip_exist = False

func_run = None
A1_arg = []
A1_bonus = []


# ========
# Parameter parsing
# ========
def extract_bonus(key, default=None, type=str, array_type=None):
    """Extract a parameter of the form "key=value" from A1_bonus, if present.

    Parameters
    ----------
    key : str
        The key part of the "key=value" bonus entry to extract

    default : ? (opt.)
        A default to return if the key is not present.

    type : type (opt.)
        Convert the extracted value to this type.

    array_type : array_type (opt.)
        Convert the extracted value to an array of this type, see array_from_string.
    """
    if key[-1] != "=":
        key += "="

    # Return the default if not present
    if not any([bonus[: len(key)] == key for bonus in A1_bonus]):
        return default

    # Extract
    for bonus in A1_bonus:
        if bonus[: len(key)] == key:
            value = bonus[len(key) :]

    # Convert
    if array_type is not None:
        return array_from_string(value, array_type)
    else:
        return type(value)


# ========
# Paths
# ========
dir_home = "/cosma/home/dp004/dc-kege1/"
dir_data5 = "/cosma5/data/dp004/dc-kege1/"
dir_data7 = "/cosma7/data/dp004/dc-kege1/"
dir_data8 = "/cosma8/data/dp217/dc-kege1/"

dir_gihr5 = dir_data5 + "gihr/"
dir_gihr7 = dir_data7 + "gihr/"
dir_gihr8 = dir_data8 + "gihr/"
dir_gihr = dir_gihr7
dir_gihrpy = dir_gihr + "gihrpy/"
dir_plots = dir_gihr + "plots/"
dir_rebound = dir_gihr + "rebound/"


# ========
# Materials
# ========
# Copy from woma for tweaking
Di_mat_id = woma.Di_mat_id
Di_mat_type = woma.Di_mat_type
type_factor = woma.type_factor

# ID offset to distinguish different bodies e.g. impactor from target
id_body = 10000
A1_mat = [mat for mat in Di_mat_id.keys()]
A1_id = [id for id in Di_mat_id.values()]
for mat, id in zip(A1_mat, A1_id):
    Di_mat_id[mat + "_2"] = id + id_body

# Default colours
Di_mat_colour = {
    # Ideal Gas
    "idg_HHe": "#1199ff",
    "idg_HHe_2": "#66ddee",
    "idg_N2": "#1199ff",
    "idg_N2_2": "#66ddee",
    "idg_CO2": "#1199ff",
    "idg_CO2_2": "#66ddee",
    # Tillotson
    "Til_iron": "#808080",
    "Til_iron_2": "#775533",
    "Til_granite": "#dd4400",
    "Til_granite_2": "#ffdd00",
    "Til_water": "#4169E1",
    "Til_water_2": "#4169E1",
    "Til_basalt": "#dd4400",
    "Til_basalt_2": "#ffdd00",
    # Hubbard & MacFarlane (1980) Uranus/Neptune
    "HM80_HHe": "#1199ff",
    "HM80_HHe_2": "#66ddee",
    "HM80_ice": "#B0C4DE",
    "HM80_ice_2": "#A080D0",
    "HM80_rock": "#708090",
    "HM80_rock_2": "#706050",
    # SESAME
    "SESAME_iron": "#808080",
    "SESAME_iron_2": "#775533",
    "SESAME_basalt": "#dd4400",
    "SESAME_basalt_2": "#ffdd00",
    "SESAME_water": "#4169E1",
    "SESAME_water_2": "#4169E1",
    "SS08_water": "#4169E1",
    "SS08_water_2": "#4169E1",
    "AQUA": "#4169E1",
    "AQUA_2": "#4169E1",
    "CMS19_H": "#1199ff",
    "CMS19_H_2": "#66ddee",
    "CMS19_He": "#1199ff",
    "CMS19_He_2": "#66ddee",
    "CD21_HHe": "#1199ff",
    "CD21_HHe_2": "#66ddee",
    # ANEOS
    "ANEOS_forsterite": "#dd4400",
    "ANEOS_forsterite_2": "#ffdd00",
    "ANEOS_iron": "#808080",
    "ANEOS_iron_2": "#775533",
    "ANEOS_Fe85Si15": "#808080",
    "ANEOS_Fe85Si15_2": "#775533",
}

# Inverted dictionaries with IDs as keys
Di_id_mat = {mat_id: mat for mat, mat_id in Di_mat_id.items()}
Di_id_type = {type_id: mat for mat, type_id in Di_mat_type.items()}
Di_id_colour = {Di_mat_id[mat]: c for mat, c in Di_mat_colour.items()}


# ========
# Particle data
# ========
# Conversion between my particle property labels and SWIFT hdf5 dataset labels
Di_hdf5_label = {  # Type
    "id": "ParticleIDs",  # L
    "mat_id": "MaterialIDs",  # i
    "pos": "Coordinates",  # d
    "vel": "Velocities",  # f
    "m": "Masses",  # f
    "h": "SmoothingLengths",  # f
    "u": "InternalEnergies",  # f
    "rho": "Densities",  # f
    "P": "Pressures",  # f
    "s": "Entropies",  # f
    "phi": "Potentials",  # f
    "fof_id": "FOFGroupIDs",  # d
}

# Friends of friends no-group ID
fof_id_none = 2**31 - 1


def load_particle_array(f, param):
    """Load a particle property array and convert units to SI.

    Parameters
    ----------
    f : h5py.File
        The opened hdf5 data file (with "r").

    param : str
        The particle property to load. See Di_hdf5_label.

    Returns
    -------
    A1_param : [?]
        The array of the particle property data (SI units).
    """
    # Number of particles
    try:
        num_picle = f["Header"].attrs["NumPart_Total"][0]
    except KeyError:
        raise Exception("Failed to load Header/NumPart_Total")

    # Units from file metadata
    file_to_SI = Conversions(
        m=float(f["Units"].attrs["Unit mass in cgs (U_M)"]) * 1e-3,
        l=float(f["Units"].attrs["Unit length in cgs (U_L)"]) * 1e-2,
        t=float(f["Units"].attrs["Unit time in cgs (U_t)"]),
    )

    # Parameter
    A1_param = np.array(f["PartType0/" + Di_hdf5_label[param]][()])

    # Convert to SI
    if param in ["id", "mat_id", "fof_id"]:
        return A1_param
    elif param == "pos":
        return (A1_param - 0.5 * f["Header"].attrs["BoxSize"]) * file_to_SI.l
    elif param == "vel":
        return A1_param * file_to_SI.v
    elif param == "m":
        return A1_param * file_to_SI.m
    elif param == "h":
        return A1_param * file_to_SI.l
    elif param in ["u", "phi"]:
        return A1_param * file_to_SI.u
    elif param == "rho":
        return A1_param * file_to_SI.rho
    elif param == "P":
        return A1_param * file_to_SI.P
    else:
        raise Exception(
            "Error: SI conversion not set for particle property: %s" % param
        )


def load_particle_data_Fp(Fp_data, A1_param):
    """Load one or more arrays of particle data directly from an hdf5 file.

    Data are returned in SI (see load_particle_array) and sorted by particle ID.

    Note: Unlike Simulation.load_snapshot_data, does not include options to
    compute derived data. So use only for data available directly from the file.

    Parameters
    ----------
    Fp_data : str
        The file path to the snapshot to load.

    A1_param : [str]
        A list of particle properties to load.

    Returns
    -------
    A1_return : [[?]]
        A list of the arrays of particle properties (SI units).
    """
    A1_return = []

    print('Loading "%s" ' % Fp_data[-52:])
    with h5py.File(Fp_data, "r") as f:
        # Load particle IDs and sort
        A1_id = load_particle_array(f, "id")
        A1_sort = np.argsort(A1_id)
        A1_id = A1_id[A1_sort]

        for param in ensure_list(A1_param):
            A1_return.append(load_particle_array(f, param)[A1_sort])

    return A1_return


# ========
# Impact scenarios
# ========
def impact_pos_vel_q_v_q_r(q, v_q, r, M_t, M_i, units_v_q="m/s", return_t=False):
    """
    Calculate the intial position and velocity of an impactor to result in the
    desired scenario at periapsis, rotated such that the velocity at periapsis
    is in the negative x direction.

    Parameters
    ----------
    q : float
        The periapsis distance (m).

    v_q : float
        The impactor's speed at periapsis, in units given by units_v_q.

    r : float
        The initial distance between body centres (m).

    M_t, M_i : float
        The masses of the target and impactor (kg).

    units_v_q : str (opt.)
        The units of the speed at periapsis: "m/s" (default), or "v_esc_q" for
        the escape speed at the periapsis distance, or "v_inf" for the speed at
        infinity instead of at periapsis (m/s), v_inf^2 = v(r)^2 - v_esc(r)^2.

    return_t : bool (opt.)
        If True, then also return the time until periapsis, e.g. to be used by
        impact_pos_vel_q_v_q_t().

    Returns
    -------
    A1_pos : [float]
        The impactor's initial x, y, z coordinates (m).

    A1_vel : [float]
        The impactor's initial x, y, z velocities (m s^-1).

    t : float (opt.)
        If return_t is True, then also return the time taken from the
        initial position until periapsis (s).
    """
    mu = G * (M_t + M_i)

    v_esc_q = np.sqrt(2 * mu / q)
    y_q = q

    # Speed at periapsis and infinity, convert if necessary
    if units_v_q == "m/s":
        v_inf = np.sqrt(v_q**2 - v_esc_q**2)
    elif units_v_q == "v_esc_q":
        v_q *= v_esc_q
        v_inf = np.sqrt(v_q**2 - v_esc_q**2)
    elif units_v_q == "v_inf":
        v_inf = v_q
        v_q = np.sqrt(v_esc_q**2 + v_inf**2)

    # Parabola
    if v_inf == 0 or v_q == v_esc_q:
        # Initial speed and position
        v = np.sqrt(2 * mu / r)
        y = v_q * y_q / v
        x = np.sqrt(r**2 - y**2)

        # True anomalies (actually the complementary angles)
        theta = np.pi - np.arccos(y_q**2 * v_q**2 / (mu * r) - 1)
    # Ellipse or hyperbola
    else:
        # Semimajor axis
        a = 1 / (2 / q - v_q**2 / mu)

        # Initial speed and position
        v = np.sqrt(mu * (2 / r - 1 / a))
        y = v_q * y_q / v
        x = np.sqrt(r**2 - y**2)

        # Eccentricity
        e = 1 - q / a

        # Check requested separation is valid
        Q = 2 * a - q
        if v_q < v_esc_q and r > Q:
            raise ValueError(
                "Invalid r = %g m for bound orbit with apoapsis = %g m" % (r, Q)
            )

        # True anomalies (actually the complementary angles)
        theta = np.arccos((1 - a * (1 - e**2) / r) / e)

    # Rotate
    if q == 0:
        phi = 0
    else:
        phi = -np.pi / 2 + theta - np.arcsin(y / r)

    x_ = x * np.cos(phi) - y * np.sin(phi)
    y_ = x * np.sin(phi) + y * np.cos(phi)
    v_x_ = -abs(v * np.cos(phi))
    v_y_ = abs(v * np.sin(phi))

    # Time until contact, if requested
    if return_t:
        # Radial
        if q == 0:
            # Parabolic
            if v_inf == 0:
                # Time until periapsis
                t = np.sqrt(2 * r**3 / (9 * mu))
            # Elliptical
            elif a > 0:
                # Standard constants
                w = 1 / r - v**2 / (2 * mu)
                wr = w * r
                # Time until periapsis
                t = (np.arcsin(np.sqrt(wr)) - np.sqrt(wr * (1 - wr))) / np.sqrt(
                    2 * mu * w**3
                )
            # Hyperbolic
            else:
                # Standard constants
                w = abs(1 / r - v**2 / (2 * mu))
                wr = w * r
                # Time until periapsis
                t = (np.sqrt(wr**2 + wr) - np.log(np.sqrt(wr) + np.sqrt(1 + wr))) / (
                    np.sqrt(2 * mu * w**3)
                )
        # Not radial
        else:
            # Parabolic
            if v_inf == 0:
                # Eccentric anomaly
                E = np.tan(0.5 * (np.pi - theta))
                # Mean anomaly
                M = E + E**3 / 3
                # Time until periapsis
                t = np.sqrt(2 * q**3 / mu) * M
            # Elliptical
            elif a > 0:
                # Eccentric anomaly
                E = np.arccos(
                    (e + np.cos(np.pi - theta)) / (1 + e * np.cos(np.pi - theta))
                )
                # Mean anomaly
                M = E - e * np.sin(E)
                # Time until periapsis
                t = np.sqrt(a**3 / mu) * M
            # Hyperbolic
            else:
                # Eccentric anomaly
                E = np.arccosh(
                    (e + np.cos(np.pi - theta)) / (1 + e * np.cos(np.pi - theta))
                )
                # Mean anomaly
                M = -E + e * np.sinh(E)
                # Time until periapsis
                t = np.sqrt(-(a**3) / mu) * M

        return np.array([x_, y_, 0]), np.array([v_x_, v_y_, 0]), t
    else:
        return np.array([x_, y_, 0]), np.array([v_x_, v_y_, 0])


def impact_pos_vel_q_v_q_t(q, v_q, t, M_t, M_i, units_v_q="m/s", r_max_factor=100):
    """
    Calculate the intial position and velocity of an impactor to result in the
    desired scenario at periapsis, rotated such that the velocity at periapsis
    is in the negative x direction.

    Find the initial distance between the body centres, r, that yields the
    desired time until contact by iteratively calling
    impact_pos_vel_q_v_q_r().

    Parameters
    ----------
    q : float
        The periapsis distance (m).

    v_q : float
        The impactor's speed at periapsis, in units given by units_v_q.

    t : float
        The time taken from the initial position until contact (s).

    M_t, M_i : float
        The masses of the target and impactor (kg).

    units_v_q : str (opt.)
        The units of the contact speed: "m/s" (default), or "v_esc_q" for the
        escape speed at the periapsis distance, or "v_inf" for the speed at
        infinity instead of at periapsis (m/s), v_inf^2 = v(r)^2 - v_esc(r)^2.

    r_max_factor : float (opt.)
        This times the sum of the body radii sets the maximum initial
        separation of the body centres for the bisection search. Default 100.

    Returns
    -------
    A1_pos : [float]
        The impactor's initial x, y, z coordinates (m).

    A1_vel : [float]
        The impactor's initial x, y, z velocities (m s^-1).
    """
    # Boundary guesses for the initial separation
    r_min = q
    r_max = r_min * r_max_factor

    # Bisection to find the separation to give the desired time to impact
    i = 0
    i_max = 1e2
    t_ = 0
    tol = 1e-6
    while tol < abs(t_ - t) / t:
        r = 0.5 * (r_min + r_max)

        try:
            t_ = impact_pos_vel_q_v_q_r(
                q, v_q, r, M_t, M_i, units_v_q=units_v_q, return_t=True
            )[2]
        except ValueError:
            # Too big r, larger than apoapsis
            t_ = 2 * t

        # Catch if r too big for the desired speed
        if np.isnan(t_):
            t_ = t * 2

        # Bisect
        if t_ < t:
            r_min = r
        else:
            r_max = r
        i += 1

        if i >= i_max:
            raise RuntimeError("Failed to find r(t) after %d iterations" % (i))

    return impact_pos_vel_q_v_q_r(q, v_q, r, M_t, M_i, units_v_q=units_v_q)


# ========
# Plotting functions
# ========
def func_cmap_param_def(po, param, A1_p):
    """Set up a colour map for a parameter.

    Parameters
    ----------
    po : PlotOptions
        An object with settings for plotting functions.

    param : str
        The parameter used to set the colours, e.g. a particle property.

    A1_p : [?]
        The parameter array (plotting units).

    po.Di_cmap_param (only required for some param)
    ----------------
    "A1_sel_unb" : [int]
        A selection of which A1_p elements correspond to unbound particles.
        Required for Q_v_inf.

    Returns
    -------
    cmap : Colormap
        The colour map.

    vmin, vmax : float
        The colour parameter limits.

    norm : Normalize
        The normalisation object.

    extend : str
        Whether to mark the bottom and/or top ends of a colour bar as extended.

    #
    # Special cases can prompt additional returns, e.g. for a second colour bar.
    #
    """
    # ========
    # Special return cases
    # ========
    if param in ["mat_id", "fof_id"]:
        # No colour map
        return None, None, None, None, None
    elif param == "Q_v_inf":
        A1_sel_unb = po.Di_cmap_param["A1_sel_unb"]

        # Q (R_unit)
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            "plasma_viridis_r",
            np.vstack(
                (
                    plt.get_cmap("plasma")(np.linspace(0.15, 0.95, 128)),
                    plt.get_cmap("viridis")(np.linspace(0.94, 0.63, 128)),
                )
            ),
        )
        vmin, vmax = po.Di_param_lim["Q"]
        norm = MidNorm(vcenter=po.R_Hill / po.r_unit, vmin=vmin, vmax=vmax)
        extend = "max"

        # v_inf (km/s)
        cmap_2 = mpl.colors.LinearSegmentedColormap.from_list(
            "", plt.get_cmap("Blues")(np.linspace(0.35, 1, 256))
        )
        if len(A1_sel_unb) > 0:
            v_inf_max = np.nanmax(A1_p[A1_sel_unb])
            if v_inf_max < 1:
                vmin_2 = 0
            else:
                vmin_2 = np.nanmin(A1_p[A1_sel_unb])
            vmax_2 = max(0.45, v_inf_max)
        else:
            vmin_2 = 0
            vmax_2 = 0.45
        norm_2 = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        extend_2 = "neither"

        return cmap, vmin, vmax, norm, extend, cmap_2, vmin_2, vmax_2, norm_2, extend_2

    # ========
    # Standard return cases
    # ========
    # Defaults
    cmap = po.cmap
    extend = "neither"
    if param in po.Di_param_lim.keys():
        vmin, vmax = po.Di_param_lim[param]
        if vmin is None:
            vmin = np.nanmin(A1_p)
        if vmax is None:
            vmax = np.nanmax(A1_p)
    else:
        vmin = np.nanmin(A1_p)
        vmax = np.nanmax(A1_p)
    norm = None

    # Specific parameter settings
    if param == "E_ext":
        cmap = cmap_plasma_viridis_r
        vmin = min(vmin, np.nanmin(A1_p))
        vmax = max(vmax, np.nanmax(A1_p))
        norm = MidNorm(vcenter=0, vmin=vmin, vmax=vmax)
        extend = "neither"
    elif param == "phase":
        cmap = ut.cmap_phase
        vmin = 0
        vmax = 6
        norm = ut.norm_phase

    # Norm, if not already set
    if norm is None:
        if param in po.Di_param_log.keys() and po.Di_param_log[param]:
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

            # Override non-positive values to the smallest positive value
            if vmin is not None and vmin <= 0:
                vmin = np.amin(A1_p[np.where(A1_p > 0)[0]])
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    return cmap, vmin, vmax, norm, extend


def func_ax_lim_most(po, time, A2_pos):
    """Set axis limits to fit around most of the particles.

    Keep a square shape.

    Parameters
    ----------
    po : PlotOptions
        An object with settings for plotting functions.

    time : float
        The snapshot time (s).

    A2_pos : [[float]]
        The particle positions (R_unit).

    Required po.Di_ax_lim
    ---------------------
    "f_extent" : float
        Fit the limits around the particles out to this fraction, e.g. 0.03.

    "pad" : float
        Pad the limits by this fraction of the axis size, e.g. 0.1.

    "pad_min" : float
        Pad the limits by at least this distance (R_unit), e.g. 0.15.

    Returns
    -------
    A1_ax_lim : [float]
        The axis limits [x_min, x_max, y_min, y_max] (R_unit).

    tick_base, tick_base_minor : [float]
        The base separation for major and minor axis ticks (R_unit).
    """
    # Find the extents of most of the particle positions
    x_min = np.quantile(A2_pos[:, 0], po.Di_ax_lim["f_extent"])
    x_max = np.quantile(A2_pos[:, 0], 1 - po.Di_ax_lim["f_extent"])
    y_min = np.quantile(A2_pos[:, 1], po.Di_ax_lim["f_extent"])
    y_max = np.quantile(A2_pos[:, 1], 1 - po.Di_ax_lim["f_extent"])

    # Make square
    dx = x_max - x_min
    dy = y_max - y_min
    x_mid = x_min + dx / 2
    y_mid = y_min + dy / 2
    if dx < dy:
        x_min = x_mid - dy / 2
        x_max = x_mid + dy / 2
    else:
        y_min = y_mid - dx / 2
        y_max = y_mid + dx / 2

    # Expand a bit further
    d = max(po.Di_ax_lim["pad"] * dx, po.Di_ax_lim["pad_min"])
    x_min -= d
    x_max += d
    y_min -= d
    y_max += d

    A1_ax_lim = np.array([x_min, x_max, y_min, y_max])

    # Nice tick spacing
    tick_base, tick_base_minor = auto_tick_spacing(A1_ax_lim[1] - A1_ax_lim[0])
    if po.tick_base_minor is None:
        tick_base_minor = None

    return A1_ax_lim, tick_base, tick_base_minor


def func_ax_lim_expand(po, time, A2_pos):
    """Set expanded axis limits to fit around roughly all of the particles.

    Stay centred on the centre of mass.

    Parameters
    ----------
    po : PlotOptions
        An object with settings for plotting functions.

    time : float
        The snapshot time (s).

    A2_pos : [[float]]
        The particle positions (R_unit).

    Required po.Di_ax_lim
    ---------------------
    "f_extent" : float
        Fit the limits around the particles out to this fraction, e.g. 5e-3.

    "pad" : float
        Pad the limits by this fraction of the axis size, e.g. 0.2.

    Returns
    -------
    A1_ax_lim : [float]
        The axis limits [x_min, x_max, y_min, y_max] (R_unit).

    tick_base, tick_base_minor : [float]
        The base separation for major and minor axis ticks (R_unit).
    """
    # Start with the input limits (assume ~centred around zero)
    A1_ax_lim = np.copy(po.A1_ax_lim).astype(float)

    # Find the extents of almost all of the particle positions
    x_min = np.quantile(A2_pos[:, 0], po.Di_ax_lim["f_extent"])
    x_max = np.quantile(A2_pos[:, 0], 1 - po.Di_ax_lim["f_extent"])
    y_min = np.quantile(A2_pos[:, 1], po.Di_ax_lim["f_extent"])
    y_max = np.quantile(A2_pos[:, 1], 1 - po.Di_ax_lim["f_extent"])

    # Centre on the centre of mass (assuming equal particle masses)
    A1_pos_com = mean(A2_pos, axis=0)
    x_com = A1_pos_com[0]
    y_com = A1_pos_com[1]

    # Scale to fit the particles
    dx_old = A1_ax_lim[1] - A1_ax_lim[0]
    dy_old = A1_ax_lim[3] - A1_ax_lim[2]
    dx_new = 2 * max(x_max - x_com, x_com - x_min, y_max - y_com, y_com - y_min)
    dx_new *= 1 + po.Di_ax_lim["pad"]

    # Enlarge axis limits (no shrinking)
    if dx_new > dx_old:
        A1_ax_lim *= dx_new / dx_old

    # Shift to center of mass
    A1_ax_lim[:2] += x_com
    A1_ax_lim[2:] += y_com

    # Placeholders
    tick_base = None
    tick_base_minor = None

    # Nice tick spacing
    tick_base, tick_base_minor = auto_tick_spacing(A1_ax_lim[1] - A1_ax_lim[0])
    if po.tick_base_minor is None:
        tick_base_minor = None

    return A1_ax_lim, tick_base, tick_base_minor


def func_ax_lim_expand_step(po, time, A2_pos):
    """Set expanded axis limits in discrete steps to fit around most of the particles.

    Parameters
    ----------
    po : PlotOptions
        An object with settings for plotting functions.

    time : float
        The snapshot time (s).

    A2_pos : [[float]]
        The particle positions (R_unit).

    Required po.Di_ax_lim
    ---------------------
    "f_extent" : float
        Fit the limits around the particles out to this fraction, e.g. 5e-3.

    "pad" : float
        Pad the limits by this fraction of the axis size, e.g. 0.2.

    "step" : float
        The step by which to increase the axis limits when needed (R_unit).

    Returns
    -------
    A1_ax_lim : [float]
        The axis limits [x_min, x_max, y_min, y_max] (R_unit).

    tick_base, tick_base_minor : [float]
        The base separation for major and minor axis ticks (R_unit).
    """
    # Start with the input limits (assume ~centred around zero)
    A1_ax_lim = np.copy(po.A1_ax_lim).astype(float)

    # Find the extents of almost all of the particle positions
    x_min = np.quantile(A2_pos[:, 0], po.Di_ax_lim["f_extent"])
    x_max = np.quantile(A2_pos[:, 0], 1 - po.Di_ax_lim["f_extent"])
    y_min = np.quantile(A2_pos[:, 1], po.Di_ax_lim["f_extent"])
    y_max = np.quantile(A2_pos[:, 1], 1 - po.Di_ax_lim["f_extent"])

    # Minimum required limits
    dx = A1_ax_lim[1] - A1_ax_lim[0]
    dy = A1_ax_lim[3] - A1_ax_lim[2]
    x_min = A1_ax_lim[1] - dx * (1 + po.Di_ax_lim["pad"])
    x_max = A1_ax_lim[0] + dx * (1 + po.Di_ax_lim["pad"])
    y_min = A1_ax_lim[3] - dy * (1 + po.Di_ax_lim["pad"])
    y_max = A1_ax_lim[2] + dy * (1 + po.Di_ax_lim["pad"])

    # Increase limits in steps to fit around the requirements
    step = po.Di_ax_lim["step"]
    while A1_ax_lim[0] > x_min:
        A1_ax_lim[0] -= step

        # To keep the same axis shape, also increase one of the y limits
        # If neither needs to be increased, then do the one closest to a required limit
        if A1_ax_lim[2] > y_min:
            A1_ax_lim[2] -= step
        elif A1_ax_lim[3] < y_min:
            A1_ax_lim[3] += step
        elif abs(A1_ax_lim[2] - y_min) < abs(A1_ax_lim[3] - y_max):
            A1_ax_lim[2] -= step
        else:
            A1_ax_lim[3] += step
    while A1_ax_lim[1] < x_max:
        A1_ax_lim[1] += step

        # As above
        if A1_ax_lim[2] > y_min:
            A1_ax_lim[2] -= step
        elif A1_ax_lim[3] < y_min:
            A1_ax_lim[3] += step
        elif abs(A1_ax_lim[2] - y_min) < abs(A1_ax_lim[3] - y_max):
            A1_ax_lim[2] -= step
        else:
            A1_ax_lim[3] += step
    while A1_ax_lim[2] > y_min:
        A1_ax_lim[2] -= step

        # As above
        if A1_ax_lim[0] > x_min:
            A1_ax_lim[0] -= step
        elif A1_ax_lim[1] < x_min:
            A1_ax_lim[1] += step
        elif abs(A1_ax_lim[0] - x_min) < abs(A1_ax_lim[1] - x_max):
            A1_ax_lim[0] -= step
        else:
            A1_ax_lim[1] += step
    while A1_ax_lim[3] < y_max:
        A1_ax_lim[3] += step

        # As above
        if A1_ax_lim[0] > x_min:
            A1_ax_lim[0] -= step
        elif A1_ax_lim[1] < x_min:
            A1_ax_lim[1] += step
        elif abs(A1_ax_lim[0] - x_min) < abs(A1_ax_lim[1] - x_max):
            A1_ax_lim[0] -= step
        else:
            A1_ax_lim[1] += step

    # Placeholders
    tick_base = None
    tick_base_minor = None

    # Nice tick spacing
    if po.tick_base is not None:
        tick_base = 0.5 * po.tick_base
        while abs(tick_base) < 0.4 * (A1_ax_lim[1] - A1_ax_lim[0]):
            tick_base *= 2

        # Keep the same relative minor-tick spacing
        if po.tick_base_minor is not None:
            tick_base_minor = tick_base * po.tick_base_minor / po.tick_base

    return A1_ax_lim, tick_base, tick_base_minor


def func_marker_size_ax_scale(po, A1_ax_lim):
    """Scale the marker size to match the axis limits.

    Parameters
    ----------
    po : PlotOptions
        An object with settings for plotting functions.

    A1_ax_lim : [float]
        The axis limits [x_min, x_max, y_min, y_max] (R_unit).

    Required po.Di_marker_size
    --------------------------
    "marker_size_ref" : float
        The reference marker size.

    "ax_size_ref" : float
        The corresponding axis size.

    "pow_size_scale" : float (opt.)
        Exponent to reduce/augment the size change with the axes change.

    Returns
    -------
    marker_size : [float]
        The scaled marker size.
    """
    ax_area = (A1_ax_lim[1] - A1_ax_lim[0]) * (A1_ax_lim[3] - A1_ax_lim[2])
    ax_area_ref = po.Di_marker_size["ax_size_ref"] ** 2

    if "pow_size_scale" in po.Di_marker_size:
        pow = po.Di_marker_size["pow_size_scale"]
    else:
        pow = 1

    return po.Di_marker_size["marker_size_ref"] * (ax_area_ref / ax_area) ** pow


# ========
# Predictions etc
# ========
def f_capture_D19(q, v_inf, M_p, R, R_Hill):
    """Estimate the captured mass fraction, following Dones (1991).

    Parameters
    ----------
    q : float
        The periapsis distance (m).

    v_inf : float
        The incoming body's speed at infinity (m/s).

    M_p : float
        The primary's mass (kg).

    R : float
        The incoming body's radius (m).
    """
    E_0 = 0.5 * v_inf**2

    DeltaE = G * M_p * R / q**2

    E_capt = -G * M_p / R_Hill

    return (0.9 * DeltaE + E_capt - E_0) / (1.8 * DeltaE)


# ========
# Misc extra functions
# ========
def plot_particle_params_init_prof(sim, param_x, param_y, ax):
    """Plot initial planet profile lines for comparison with particle parameters."""
    po = sim.po

    if sim.impact is not None:
        ip_1 = sim.impact.ic_1.ip
        ip_2 = sim.impact.ic_2.ip
        ip_1.load_profiles()
        ip_2.load_profiles()
        A1_ip = [ip_1, ip_2]
    else:
        ip = sim.ic.ip
        ip.load_profiles()
        A1_ip = [ip]

    colour = "k"
    alpha = 0.9
    lw = 0.7

    # Plot the profile for each body and material
    for ip_i in A1_ip:
        for mat_id in ip_i.A1_mat_id:
            A1_sel_mat_ip_i = np.where(ip_i.A1_mat_id == mat_id)[0]
            if param_x in ["r", "r_xy"]:
                A1_x_prof = ip_i.A1_r[A1_sel_mat_ip_i] / po.r_unit
            elif param_x == "rho":
                A1_x_prof = ip_i.A1_rho[A1_sel_mat_ip_i]
            elif param_x == "T":
                A1_x_prof = ip_i.A1_T[A1_sel_mat_ip_i]
            elif param_x == "P":
                A1_x_prof = ip_i.A1_P[A1_sel_mat_ip_i]
            elif param_x == "u":
                A1_x_prof = ip_i.A1_u[A1_sel_mat_ip_i]
            elif param_x == "s":
                A1_x_prof = woma.A1_s_rho_T(
                    ip_i.A1_rho[A1_sel_mat_ip_i],
                    ip_i.A1_T[A1_sel_mat_ip_i],
                    ip_i.A1_mat_id[A1_sel_mat_ip_i],
                )

            if param_y in ["r", "r_xy"]:
                A1_y_prof = ip_i.A1_r[A1_sel_mat_ip_i] / po.r_unit
            elif param_y == "rho":
                A1_y_prof = ip_i.A1_rho[A1_sel_mat_ip_i]
            elif param_y == "T":
                A1_y_prof = ip_i.A1_T[A1_sel_mat_ip_i]
            elif param_y == "P":
                A1_y_prof = ip_i.A1_P[A1_sel_mat_ip_i]
            elif param_y == "u":
                A1_y_prof = ip_i.A1_u[A1_sel_mat_ip_i]
            elif param_y == "s":
                A1_y_prof = woma.A1_s_rho_T(
                    ip_i.A1_rho[A1_sel_mat_ip_i],
                    ip_i.A1_T[A1_sel_mat_ip_i],
                    ip_i.A1_mat_id[A1_sel_mat_ip_i],
                )

            plt.plot(A1_x_prof, A1_y_prof, c=colour, lw=lw, alpha=alpha)


def plot_particle_params_init_orbits(sim, param_x, param_y, ax):
    """Plot the parameters of initial orbits for comparison with particle data."""
    o_0_1, o_0_2, o_0_c = sim.init_orbits()

    # Plot each pre-impact body and the centre of mass
    for i_o, (o_i, c_i) in enumerate([[o_0_1, c_1], [o_0_2, c_2], [o_0_c, c_c]]):
        if param_x == "a":
            x_i = o_i.a / po.r_unit
        elif param_x == "e":
            x_i = o_i.e
        elif param_x == "q":
            x_i = o_i.q / po.r_unit
        elif param_x == "Q":
            x_i = o_i.Q / po.r_unit

        if param_y == "a":
            y_i = o_i.a / po.r_unit
        elif param_y == "e":
            y_i = o_i.e
        elif param_y == "i":
            y_i = o_i.i * rad_to_deg * np.cos(o_i.Omega - pi / 2)
        elif param_y == "q":
            y_i = o_i.q / po.r_unit
        elif param_y == "Q":
            y_i = o_i.Q / po.r_unit
        elif param_y == "a_eq":
            y_i = o_i.a * (1 - o_i.e**2) * np.cos(o_i.i) ** 2 / po.r_unit

        ax.scatter(
            [x_i],
            [y_i],
            marker="o",
            facecolor=c_i,
            edgecolor="k",
            s=8**2,
            lw=0.5,
            alpha=0.8,
            zorder=9,
        )


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR utilities.py  ====\n")
