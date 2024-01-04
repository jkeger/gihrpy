"""
More classes for organising GIHR projects, specifically for main simulations.

See classes.py.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.utilities_phodei as ut_ph
from gihrpy.classes import PlotOptions, InitProf, InitCond


class ImpactInitCond:
    """Initial conditions for an impact scenario (or similar).

    Parameters
    ----------
    copy : str
        A name of a Simulation object that has an ImpactInitCond object from
        which to copy all attributes that are not provided as arguments here.

    ic_1, ic_2 : InitCond
        The initial conditions objects for each of the two bodies.

    Fp_init_cond_1, Fp_init_cond_2 : str or (Sim, int)
        The file path for the saved initial conditions for each body, or a tuple
        of the name of another Simulation object and the snapshot time.

    A1_preset_1, A1_preset_2 : [str]
        Options for tweaks to the initial particle sets:
        "v=0"           Set all particle velocities to zero (unless spinning).
        "rotate=AB"     Rotate the particles about the A=x,y,z axis by B degrees.
        "spin=A"        Shorthand rotate for spinning bodies so the angular
                        momentum is in the A-axis direction.
        "xp"            Use an external potential instead, so no particles.
        "fof=A"         Extract the fof group with ID=A instead of all particles.

    A1_pos_1_in, A1_vel_1_in, A1_pos_2_in, A1_vel_2_in : [float]
        The input cartesian positions and velocities of the two bodies, before
        any reference-frame changes (m, m s^-1).

    b, B, q : float (only one)
    v_c, v_c_esc, v_inf, v_q, e : float (only one)
    t_c, t_q : float (only one)
        Parameters to override the input positions and velocities, to pass to
        woma.impact_pos_vel_b_v_c_t() or ut.impact_pos_vel_q_v_q_t(),
        to calculate A1_pos_2_in and A1_vel_2_in with A1_pos_1_in and
        A1_vel_1_in set to all zeros.
        b           Impact parameter.
        B           Impact angle (rad).
        q           Periapsis (m).
        v_c         Speed at contact (m s^-1).
        v_c_esc     Speed at contact (v_esc).
        v_inf       Speed at infinity (m s^-1).
        v_q         Speed at periapsis (m s^-1).
        e           Eccentricity.
        t_c         Time to contact (s).
        t_q         Time to periapsis (s).

    is_centre_mass : bool
        If True, then shift coordinates to the centre-of-mass frame.

    is_centre_mom : bool
        If True, then shift velocities to the zero-momentum frame.

    m_xp : float
        The mass of an external potential point mass (kg).

    A1_pos_xp : [float]
        The input cartesian positions of an external potential point mass,
        before any reference-frame changes (m).

    A1_pos_com, A1_vel_com : [float], [float]
        The input cartesian positions and velocities of the centre of mass
        of the particles, e.g. around an external potential point mass,
        before any reference-frame changes (m, m s^-1).

    time_c : float
        The approximate time of a collision from a rebound simulation (yr), with
        which most other parameters will be overridden.

    Attributes
    ----------
    M_1, M_2, M_tot : float
        The masses of the two bodies and the total mass (kg).

    A1_pos_1, A1_vel_1, A1_pos_2, A1_vel_2 : [float]
        The final cartesian positions and velocities of the two bodies, after
        any reference-frame changes (m, m s^-1).

    num_picle : int
        The total number of particles.

    m_picle : float
        The average particle mass (kg).

    v_esc : float
        The two-body mutual escape speed (m s^-1).

    A1_L : float
        The orbital angular momentum of the system (km m^2 s^-1).
    """

    def __init__(
        self,
        copy=None,
        ic_1=None,
        ic_2=None,
        Fp_init_cond_1=None,
        Fp_init_cond_2=None,
        A1_preset_1=None,
        A1_preset_2=None,
        A1_pos_1_in=None,
        A1_vel_1_in=None,
        A1_pos_2_in=None,
        A1_vel_2_in=None,
        b=None,
        B=None,
        q=None,
        v_c=None,
        v_c_esc=None,
        v_inf=None,
        v_q=None,
        e=None,
        t_c=None,
        t_q=None,
        is_centre_mass=None,
        is_centre_mom=None,
        m_xp=None,
        A1_pos_xp=None,
        A1_pos_com=None,
        A1_vel_com=None,
        time_c=None,
    ):
        self.copy = copy
        self.ic_1 = ic_1
        self.ic_2 = ic_2
        self.Fp_init_cond_1 = Fp_init_cond_1
        self.Fp_init_cond_2 = Fp_init_cond_2
        self.A1_preset_1 = A1_preset_1
        self.A1_preset_2 = A1_preset_2
        self.A1_pos_1_in = A1_pos_1_in
        self.A1_vel_1_in = A1_vel_1_in
        self.A1_pos_2_in = A1_pos_2_in
        self.A1_vel_2_in = A1_vel_2_in
        self.b = b
        self.B = B
        self.q = q
        self.v_c = v_c
        self.v_c_esc = v_c_esc
        self.v_inf = v_inf
        self.v_q = v_q
        self.e = e
        self.t_c = t_c
        self.t_q = t_q
        self.is_centre_mass = is_centre_mass
        self.is_centre_mom = is_centre_mom
        self.m_xp = m_xp
        self.A1_pos_xp = A1_pos_xp
        self.A1_pos_com = A1_pos_com
        self.A1_vel_com = A1_vel_com
        self.time_c = time_c

        self.init_done = False

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        # Input snapshot files
        if type(self.Fp_init_cond_1) == tuple:
            sim_1, time_1 = self.Fp_init_cond_1
            sim_1.init()
            self.Fp_init_cond_1 = sim_1.A1_Fp_snap[time_1]
        if type(self.Fp_init_cond_2) == tuple:
            sim_2, time_2 = self.Fp_init_cond_2
            sim_2.init()
            self.Fp_init_cond_2 = sim_2.A1_Fp_snap[time_2]

        # Set converted values to ensure no copy overwrites
        if self.b is not None and self.B is None:
            self.B = np.arcsin(self.b)
        elif self.B is not None and self.b is None:
            self.b = np.sin(self.B)

        # Exclude copying velocities before converting, in case e.g. need to copy v_esc
        A1_exclude = []
        if self.t_c is not None and self.v_c is None and self.v_c_esc is None:
            A1_exclude.append("v_c")
            A1_exclude.append("v_c_esc")
        if self.v_c is not None and self.v_c_esc is None:
            A1_exclude.append("v_c_esc")
        elif self.v_c_esc is not None and self.v_c is None:
            A1_exclude.append("v_c")

        # Copy
        if self.copy is not None:
            sim = Simulation._Di_simulation[self.copy]
            sim.init()
            copy = sim.impact
            copy.init()
            copy_object(self, copy, A1_set_eq=["ic_1", "ic_2"], A1_exclude=A1_exclude)

        self.ic_1.init()
        self.ic_2.init()

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        # Skip the rest if derived from a rebound simulation
        if self.time_c is not None:
            # Placeholders
            self.num_picle_1 = -1
            self.num_picle_2 = -1
            self.num_picle = -1
            return

        # Particle numbers and masses
        if self.A1_preset_1 is not None and "xp" in self.A1_preset_1:
            self.M_1 = self.m_xp
            self.num_picle_1 = 0
        else:
            self.M_1 = self.ic_1.M
            self.num_picle_1 = self.ic_1.num_picle
        self.num_picle_2 = self.ic_2.num_picle
        self.num_picle = self.num_picle_1 + self.num_picle_2
        self.M_2 = self.ic_2.M
        self.M_tot = self.M_1 + self.M_2
        if self.m_xp is not None:
            self.m_picle = self.M_2 / self.num_picle
        else:
            self.m_picle = self.M_tot / self.num_picle

        # Speed units
        self.v_esc = v_esc_from_M_R(self.M_1 + self.M_2, self.ic_1.R_s + self.ic_2.R_s)
        if self.t_c is not None and self.v_c is None and self.v_c_esc is None:
            # Speed at contact from v_inf
            r_c = self.ic_1.R_s + self.ic_2.R_s
            self.v_c = np.sqrt(self.v_inf**2 + 2 * G * self.M_tot / r_c)
        if self.v_c is not None and self.v_c_esc is None:
            self.v_c_esc = self.v_c / self.v_esc
        elif self.v_c_esc is not None and self.v_c is None:
            self.v_c = self.v_c_esc * self.v_esc

        # Ensure array types
        if self.A1_pos_xp is not None:
            self.A1_pos_xp = np.array(self.A1_pos_xp)

    def prep_scenario(self):
        """Compute and set the initial positions and velocities, etc."""
        # Compute input positions and velocities
        if self.t_c is not None or self.t_q is not None:
            # Use the spherical surface radii for spinning planets to keep the
            # same scenario regardless of spin
            if self.ic_1.ip.period is None:
                R_s_1 = self.ic_1.R_s
            else:
                R_s_1 = self.ic_1.R_s_sph
            if self.ic_2.ip.period is None:
                R_s_2 = self.ic_2.R_s
            else:
                R_s_2 = self.ic_2.R_s_sph

            self.A1_pos_1_in, self.A1_vel_1_in = [0, 0, 0], [0, 0, 0]
            # By angle and speed at contact or infinity
            if self.b is not None:
                if self.v_c is not None:
                    self.A1_pos_2_in, self.A1_vel_2_in = woma.impact_pos_vel_b_v_c_t(
                        self.b,
                        self.v_c,
                        self.t_c,
                        R_s_1,
                        R_s_2,
                        self.ic_1.M,
                        self.ic_2.M,
                        units_b="b",
                        units_v_c="m/s",
                    )
                elif self.v_inf is not None:
                    self.A1_pos_2_in, self.A1_vel_2_in = woma.impact_pos_vel_b_v_c_t(
                        self.b,
                        self.v_inf,
                        self.t_c,
                        R_s_1,
                        R_s_2,
                        self.ic_1.M,
                        self.ic_2.M,
                        units_b="b",
                        units_v_c="v_inf",
                    )
            # By periapsis and speed at periapsis or infinity
            elif self.q is not None:
                # First get periapsis speed from eccentricity
                if self.e is not None:
                    mu = G * (self.M_1 + self.M_2)
                    if self.e == 1:
                        self.v_q = np.sqrt(2 * mu / self.q)
                    else:
                        a = self.q / (1 - self.e)
                        self.v_q = np.sqrt(mu * (2 / self.q - 1 / a))

                if self.v_q is not None:
                    self.A1_pos_2_in, self.A1_vel_2_in = ut.impact_pos_vel_q_v_q_t(
                        self.q,
                        self.v_q,
                        self.t_q,
                        self.M_1,
                        self.M_2,
                        units_v_q="m/s",
                    )
                elif self.v_inf is not None:
                    self.A1_pos_2_in, self.A1_vel_2_in = ut.impact_pos_vel_q_v_q_t(
                        self.q,
                        self.v_inf,
                        self.t_q,
                        self.M_1,
                        self.M_2,
                        units_v_q="v_inf",
                    )
                else:
                    raise Exception
            else:
                raise Exception

        # Ensure array types
        self.A1_pos_1_in = np.array(self.A1_pos_1_in, dtype=float)
        self.A1_vel_1_in = np.array(self.A1_vel_1_in, dtype=float)
        self.A1_pos_2_in = np.array(self.A1_pos_2_in, dtype=float)
        self.A1_vel_2_in = np.array(self.A1_vel_2_in, dtype=float)

        # Positions and velocities
        self.A1_pos_1 = self.A1_pos_1_in.copy()
        self.A1_vel_1 = self.A1_vel_1_in.copy()
        self.A1_pos_2 = self.A1_pos_2_in.copy()
        self.A1_vel_2 = self.A1_vel_2_in.copy()

        self.A1_pos_com = centre_of_mass(
            [self.M_1, self.M_2], [self.A1_pos_1, self.A1_pos_2]
        )
        self.A1_vel_com = centre_of_mass(
            [self.M_1, self.M_2], [self.A1_vel_1, self.A1_vel_2]
        )

        # Orbital angular momentum in CoM frame
        self.A1_L = self.M_1 * np.cross(
            self.A1_pos_1 - self.A1_pos_com, self.A1_vel_1 - self.A1_vel_com
        ) + self.M_2 * np.cross(
            self.A1_pos_2 - self.A1_pos_com, self.A1_vel_2 - self.A1_vel_com
        )

        # Change to centre of mass and/or momentum
        if self.is_centre_mass:
            self.A1_pos_1 -= self.A1_pos_com
            self.A1_pos_2 -= self.A1_pos_com

        if self.is_centre_mom:
            self.A1_vel_1 -= self.A1_vel_com
            self.A1_vel_2 -= self.A1_vel_com


class Simulation:
    """All relevant information for a SWIFT impact or other simulation.

    Parameters
    ----------
    name : str
        The simulation name.

    copy : str
        A name of a different Simulation object from which to copy all
        attributes that are not provided as other arguments here.

    ic : InitCond
        An object with informantion about the initial conditions, for a settling
        simulation.

    impact : ImpactInitCond
        An object with information about the initial conditions, for an impact
        simulation.

    po : PlotOptions
        An object with settings for plotting functions. Defaults to inherit from
        ic if not provided.

    dir_proj : str
        The path to the base directory for the project, expected to contain a
        sub directory for this simulation. Defaults to copy ic or impact.ic_1.

    category : str
        The base type of simulation, e.g. the project name.

    A1_time : [int]
        A list of all snapshot times (s).

    A1_time_snap, A1_time_snip : [int]
        Alternatively to A1_time, the list of snapshot times excluding snipshot
        times (s), and vice versa. Any values in A1_time_snip that also appear
        in A1_time_snap are removed from A1_time_snip. A1_time is then set to be
        the concatenation of the two. If not provided, then A1_time_snap
        defaults to A1_time.

    file_to_SI : Conversions
        The unit conversion object for the simulation code's units into SI
        units. See jkeger.py.

    boxsize : float
        The extent of the simulation box (m).

    h_max : float
        The maximum smoothing length (m).

    soft : float
        The gravitational softening (m), typically roughly equal to the inter-
        particle separation.

    link_len : str
        The friends of friends linking length (file_to_SI.l), as a string to
        ensure a match to the SWIFT output filenames.

    A1_link_len : [str]
        A list of additional linking lengths, e.g. one lower and one higher for
        constraining the sensitivity.

    snap_id_type : str
        "time"      Time labels set by an output list (default).
        "sequence"  SWIFT's default 4-digit padded count IDs.

    Fp_init_cond : str
        The file path for the initial conditions. Defaults to from ic if
        provided, or for an impact to <dir_proj>/init_cond/init_<name>.hdf5.

    rsim_in : RebSim
        A rebound simulation object to use for the initial conditions.

    A1_copy_extra : [str]
        If provided, then set or overwrite the listed attributes from the copy
        object after normal setup, e.g. Fp_init_cond or derived attributes.

    Di_po_edit : {}
        A dictionary of parameter names and values to edit the PlotOptions
        object.

    Di_misc : {}
        A dictionary of miscellaneous parameters.

    Attributes
    ----------
    num_picle : int
        The total number of particles.

    A1_snap_id : [str]
        The list of snapshot IDs as string integers with padded zeros.

    dir_snap : str
        The path to the directory for snapshots and other data files:
        dir_proj/name/snapshots/.

    Fp_snap_stem : str
        The stem file path for snapshots and other data files:
        dir_proj/name/snapshots/name.

    A1_Fp_snap : [str]
        The list of full snapshot file paths.

    Fp_accum_data : str
        The file path for accumulated snapshot etc data.

    A1_mat_sim, A1_mat_id_sim : [str], [int]
        The material names and IDs present in this simulation.
    """

    # Dictionary of all Simulation objects
    _Di_simulation = {}

    def __init__(
        self,
        name,
        copy=None,
        ic=None,
        impact=None,
        po=None,
        dir_proj=None,
        category=None,
        A1_time=None,
        A1_time_snap=None,
        A1_time_snip=None,
        file_to_SI=None,
        boxsize=None,
        h_max=None,
        soft=None,
        link_len=None,
        A1_link_len=None,
        snap_id_type=None,
        Fp_init_cond=None,
        rsim_in=None,
        A1_copy_extra=None,
        Di_po_edit=None,
        Di_misc=None,
    ):
        self.name = name
        self.copy = copy
        self.ic = ic
        self.impact = impact
        self.po = po
        self.dir_proj = dir_proj
        self.category = category
        self.A1_time = A1_time
        self.A1_time_snap = A1_time_snap
        self.A1_time_snip = A1_time_snip
        self.file_to_SI = file_to_SI
        self.boxsize = boxsize
        self.h_max = h_max
        self.soft = soft
        self.link_len = link_len
        self.A1_link_len = A1_link_len
        self.snap_id_type = snap_id_type
        self.Fp_init_cond = Fp_init_cond
        self.rsim_in = rsim_in
        self.A1_copy_extra = A1_copy_extra
        self.Di_po_edit = Di_po_edit
        self.Di_misc = Di_misc

        self.init_done = False

        # Add to the dictionary
        if self.name is not None:
            Simulation._Di_simulation[self.name] = self

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        if self.ic is not None:
            self.ic.init()
        if self.impact is not None:
            self.impact.init()

        # Overrides (before copy)
        if self.A1_time is None and self.A1_time_snap is not None:
            # Derive all times from snapshot and snipshot times
            assert self.A1_time_snip is not None
            # Remove duplicates
            self.A1_time_snip = np.setdiff1d(self.A1_time_snip, self.A1_time_snap)
            # All times
            self.A1_time = np.unique(
                np.concatenate((self.A1_time_snap, self.A1_time_snip))
            )
        elif self.A1_time_snap is None and self.A1_time is not None:
            # Derive snapshot times from all times, assuming no snipshots
            self.A1_time_snap = self.A1_time
            self.A1_time_snip = np.array([])

        # Reset defaults to None to not affect copying
        if self.Di_misc == {}:
            self.Di_misc = None

        # Copy
        if self.copy is not None:
            copy = Simulation._Di_simulation[self.copy]
            copy.init()
            copy_object(
                self,
                copy,
                A1_set_eq=["ic", "impact", "po"],
                A1_exclude=["Fp_init_cond", "A1_copy_extra"],
            )

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        if self.rsim_in is not None:
            self.rsim_in.init()

        # Defaults
        if self.po is None:
            if self.ic is not None:
                self.po = self.ic.po
            else:
                self.po = PlotOptions()
        if self.Di_misc is None:
            self.Di_misc = {}

        # Derived attributes
        self.SI_to_file = self.file_to_SI.inv()

        # File paths
        if self.dir_proj is None:
            if self.ic is not None:
                self.dir_proj = self.ic.dir_proj
            elif self.impact is not None:
                self.dir_proj = self.impact.ic_1.dir_proj
        self.dir_snap = "%s/%s/snapshots/" % (self.dir_proj, self.name)
        self.Fp_snap_stem = "%s/%s" % (self.dir_snap, self.name)

        if self.Fp_init_cond is None:
            if self.ic is not None:
                self.Fp_init_cond = self.ic.Fp_save
            elif self.impact is not None:
                self.Fp_init_cond = "%s/init_cond/init_%s.hdf5" % (
                    self.dir_proj,
                    self.name,
                )

        if self.snap_id_type == "sequence":
            self.A1_snap_id = np.array(["%04d" % i for i in range(len(self.A1_time))])
        else:
            self.A1_snap_id = np.array(["%d" % time for time in self.A1_time])
        self.A1_Fp_snap = np.array(
            ["%s_%s.hdf5" % (self.Fp_snap_stem, snap_id) for snap_id in self.A1_snap_id]
        )
        self.Fp_accum_data = "%s_accum_data.hdf5" % (self.Fp_snap_stem)

        # Include initial conditions in snapshot list
        if self.A1_time[0] != 0:
            self.A1_time = np.append(0, self.A1_time)
            self.A1_time_snap = np.append(0, self.A1_time_snap)
            if self.snap_id_type == "sequence":
                self.A1_snap_id = np.append("0000", self.A1_snap_id)
            else:
                self.A1_snap_id = np.append("0", self.A1_snap_id)
            self.A1_Fp_snap = np.append(self.Fp_init_cond, self.A1_Fp_snap)
        elif self.impact is not None or len(self.A1_time) == 1:
            self.A1_Fp_snap[0] = self.Fp_init_cond

        # Number(s) of particles
        if self.impact is not None:
            # Set extra info if derived from a rebound simulation
            if self.Di_misc is not None and "N_c_1" in self.Di_misc.keys():
                self.impact.num_picle_1 = self.Di_misc["N_c_1"]
                self.impact.num_picle_2 = self.Di_misc["N_c_2"]
                self.impact.num_picle = (
                    self.impact.num_picle_1 + self.impact.num_picle_2
                )
        if self.ic is not None:
            self.num_picle = self.ic.num_picle
        elif self.impact is not None:
            self.num_picle = self.impact.num_picle

        # Materials
        if self.impact is not None:
            self.A1_mat_id_sim = np.append(
                self.impact.ic_1.ip.A1_mat_id_layer,
                np.array(self.impact.ic_2.ip.A1_mat_id_layer) + ut.id_body,
            )
        else:
            self.A1_mat_id_sim = np.array(self.ic.ip.A1_mat_id_layer)
        self.A1_mat_sim = [ut.Di_id_mat[mat_id] for mat_id in self.A1_mat_id_sim]
        if self.A1_link_len == "None":
            self.A1_link_len = None

        # Overrides etc
        if self.po.ref_m_picle is not None:
            # Scale the marker size according to the particle mass, e.g. for
            # consistency across resolutions
            m_rel = self.m_picle / self.po.ref_m_picle
            if self.po.marker_size is not None:
                self.po.marker_size = (
                    np.sqrt(self.po.marker_size) * np.cbrt(m_rel)
                ) ** 2
        if self.A1_copy_extra is not None:
            for var in self.A1_copy_extra:
                exec("self.%s = copy.%s" % (var, var))
        if self.Di_po_edit is not None:
            # Make a clean copy to edit
            self.po = deepcopy(self.po)

            for var in self.Di_po_edit.keys():
                # self.po.var = self.Di_po_edit["var"]
                exec("self.po." + var + ' = self.Di_po_edit["' + var + '"]')

        # Placeholders
        self.A1_id_0 = None
        self.A1_m_0 = None
        self.A1_mat_id_0 = None

    # ========
    # Data loading etc
    # ========
    def load_or_compute_oset(self, time, A2_pos, A2_vel, M_p, A1_m, do_recalc=False):
        """Compute and save particle orbits, or load if previously saved.

        Parameters
        ----------
        time : int
            The snapshot time (s).

        A2_pos, A2_vel : [float]
            The particle positions (m) and velocities (m/s) relative to the
            primary.

        M_p : float
            The primary mass (kg).

        A1_m : [float]
            The particle masses (kg).

        do_recalc : bool (opt.)
            Force recalculation even if a previously saved file exists.

        Returns
        -------
        oset : OrbitSet
            The set of particle orbits.
        """
        self.init()
        Fp_oset = "%s_oset_%06d.pkl" % (self.Fp_snap_stem, time)

        # Do (re)calculation if the file doesn't exist or was made before an update
        if not os.path.isfile(Fp_oset):
            do_recalc = True
        else:
            time_mod = os.path.getmtime(Fp_oset)
            time_upd = mktime(
                datetime.datetime.strptime("2022/12/01", "%Y/%m/%d").timetuple()
            )
            if time_mod < time_upd:
                do_recalc = True

        # Compute and save, or load
        if do_recalc:
            print(
                "Computing %s %d orbits... " % (self.name[-48:], time),
                flush=True,
                end="",
            )
            oset = OrbitSet(
                A2_pos=A2_pos,
                A2_vel=A2_vel,
                A1_M_p=np.full(len(A1_m), M_p),
                A1_m=A1_m,
            )
            print("Done")

            print('\rWriting "%s" ' % Fp_oset[-52:], flush=True, end="")
            with open(Fp_oset, "wb") as f_oset:
                pickle.dump(oset, f_oset)
            print("Done      ")
        else:
            print('Loading "%s" ' % Fp_oset[-58:], flush=True, end="")
            with open(Fp_oset, "rb") as f_oset:
                oset = pickle.load(f_oset)
            print("Done")

        return oset

    def load_or_compute_oset_fof(
        self, time, link_len=None, fof_id_max=None, do_recalc=False
    ):
        """Compute and save fof-group orbits, or load if previously saved.

        Parameters
        ----------
        time : int
            The snapshot time (s).

        link_len : str (opt.)
            The fof linking length (file_to_SI.l).

        fof_id_max : int (opt.)
            Only include fof groups up to and including this ID.

        do_recalc : bool (opt.)
            Force recalculation even if a previously saved file exists.

        Returns
        -------
        oset_fof : OrbitSet
            The set of fof-group orbits.
        """
        self.init()

        # Defaults
        if link_len is None:
            link_len = self.link_len
        if fof_id_max is None:
            fof_id_max = self.po.fof_id_max

        Fp_oset_fof = "%s_oset_fof_%06d_%s.pkl" % (self.Fp_snap_stem, time, link_len)

        # Do (re)calculation if the file doesn't exist or was made before an update
        if not os.path.isfile(Fp_oset_fof):
            do_recalc = True
        else:
            time_mod = os.path.getmtime(Fp_oset_fof)
            time_upd = mktime(
                datetime.datetime.strptime("2022/12/01", "%Y/%m/%d").timetuple()
            )
            if time_mod < time_upd:
                do_recalc = True

        # Compute and save, or load
        if do_recalc:
            # Load snapshot data
            A1_m, A2_pos, A2_vel, A1_fof_id = self.load_snapshot_data(
                time, ["m", "pos", "vel", "fof_id"], link_len=link_len
            )

            print(
                "Computing %s %d fof orbits... " % (self.name[-48:], time),
                flush=True,
                end="",
            )

            # Compute each group's orbit
            A1_o = []
            for fof_id in range(fof_id_max + 1):
                # Select the particles in the group
                A1_sel_fof = np.where(A1_fof_id == fof_id)[0]

                # Group centre of mass position and velocity
                A1_pos = centre_of_mass(A1_m[A1_sel_fof], A2_pos[A1_sel_fof])
                A1_vel = centre_of_mass(A1_m[A1_sel_fof], A2_vel[A1_sel_fof])
                m = sum(A1_m[A1_sel_fof])

                # Estimate the orbital parameters
                A1_o.append(
                    Orbit(
                        A1_pos - self.A1_pos_p,
                        A1_vel - self.A1_vel_p,
                        self.M_p,
                        m,
                        id=fof_id,
                    )
                )

            oset_fof = OrbitSet(A1_o=A1_o)
            print("Done")

            print('\rWriting "%s" ' % Fp_oset_fof[-58:], flush=True, end="")
            with open(Fp_oset_fof, "wb") as f_oset_fof:
                pickle.dump(oset_fof, f_oset_fof)
            print("Done")
        else:
            print('Loading "%s" ' % Fp_oset_fof[-58:], flush=True, end="")
            with open(Fp_oset_fof, "rb") as f_oset_fof:
                oset_fof = pickle.load(f_oset_fof)

            # Force recalculation if the required groups are not in the saved set
            if fof_id_max != ut.fof_id_none and np.amax(oset_fof.A1_id) < fof_id_max:
                self.load_or_compute_oset_fof(
                    time, link_len=link_len, fof_id_max=fof_id_max, do_recalc=True
                )

            print("Done")

        return oset_fof

    def load_snapshot_data(self, time, A1_param, link_len=None):
        """Load (or compute) one or more arrays of snapshot particle data.

        Parameters
        ----------
        time : int
            The time of the snapshot to load (s).
            0       The initial conditions.
            -1      The final snapshot.

        A1_param : [str]
            A list of particle properties to load, or one. Non-standard options:
            o           The OrbitSet of particle orbits.
            Q_v_inf     Apoapsis if bound, speed at infinity if unbound.

        link_len : str
            The fof linking length (file_to_SI.l), required to load fof_id.

        Returns
        -------
        A1_return : [[?]]
            A list of the arrays of particle properties, or one.
        """
        self.init()
        time = int(time)

        A1_return = []
        A1_param = ensure_list(A1_param)

        Fp_snap = self.Fp_snap_from_time(time)

        print('Loading "%s" ' % Fp_snap[-64:])
        with h5py.File(Fp_snap, "r") as f:
            # Load particle IDs and sort
            A1_id = ut.load_particle_array(f, "id")
            A1_sort = np.argsort(A1_id)
            A1_id = A1_id[A1_sort]

            # Parameters with shared requirements
            A1_param_pos = ["pos", "r", "x", "y", "z"]
            A1_param_vel = ["vel", "v", "v_x", "v_y", "v_z"]
            A1_param_orbit = ["o", "a", "e", "i", "q", "Q", "a_eq", "Q_v_inf"]

            # ========
            # Load (and sort) any auxiliary data needed to compute the requested data
            # ========
            # Initial conditions masses and material IDs
            if any(
                param in ["m", "mat_id", "T", "s"] + A1_param_orbit
                for param in A1_param
            ):
                # Load if not already loaded
                if self.A1_id_0 is None:
                    print('Loading "%s" ' % self.Fp_init_cond[-64:])
                    with h5py.File(self.Fp_init_cond, "r") as f_0:
                        A1_id_0 = ut.load_particle_array(f_0, "id")
                        A1_sort_0 = np.argsort(A1_id_0)
                        A1_id_0 = A1_id_0[A1_sort_0]

                        A1_m_0 = ut.load_particle_array(f_0, "m")[A1_sort_0]
                        A1_mat_id_0 = ut.load_particle_array(f_0, "mat_id")[A1_sort_0]

                    # Store for future reference
                    self.A1_id_0 = A1_id_0
                    self.A1_m_0 = A1_m_0
                    self.A1_mat_id_0 = A1_mat_id_0

                # Select the initial particles that match the particles in this snapshot
                A1_match = np.intersect1d(A1_id, self.A1_id_0, assume_unique=True)
                A1_sel = np.searchsorted(self.A1_id_0, A1_match)

            # Masses
            if any(param in ["m"] + A1_param_orbit for param in A1_param):
                A1_m = deepcopy(self.A1_m_0)[A1_sel]

            # Material IDs
            if any(param in ["mat_id", "T", "s"] for param in A1_param):
                A1_mat_id = deepcopy(self.A1_mat_id_0)[A1_sel]

                # Identify the particles in the second body
                if self.impact is not None:
                    A1_sel_2 = np.where(self.impact.num_picle_1 <= A1_id)[0]
                    A1_mat_id[A1_sel_2] += ut.id_body

            # Positions
            if any(param in A1_param_pos + A1_param_orbit for param in A1_param):
                A2_pos = ut.load_particle_array(f, "pos")[A1_sort]

            # Velocities
            if any(param in A1_param_vel + A1_param_orbit for param in A1_param):
                A2_vel = ut.load_particle_array(f, "vel")[A1_sort]

            # Densities and internal energies
            if any(param in ["rho", "u", "T", "s"] for param in A1_param):
                A1_rho = ut.load_particle_array(f, "rho")[A1_sort]
                A1_u = ut.load_particle_array(f, "u")[A1_sort]

                # Prep EoS tables for each material
                if any(param in ["T", "s"] for param in A1_param):
                    for mat_id in np.unique(A1_mat_id % ut.id_body):
                        woma.load_eos_tables(Di_id_mat[mat_id])

            # Orbits
            if any(param in A1_param_orbit for param in A1_param):
                oset = self.load_or_compute_oset(
                    time, A2_pos - self.A1_pos_p, A2_vel - self.A1_vel_p, self.M_p, A1_m
                )

            # ========
            # Load or compute each requested array
            # ========
            for param in A1_param:
                if param == "id":
                    A1_data = A1_id

                elif param == "m":
                    A1_data = A1_m

                elif param == "mat_id":
                    A1_data = A1_mat_id

                elif param == "pos":
                    A1_data = A2_pos

                elif param == "r":
                    A1_data = root_sum_sq(A2_pos.T)

                elif param == "x":
                    A1_data = A2_pos[:, 0]

                elif param == "y":
                    A1_data = A2_pos[:, 1]

                elif param == "z":
                    A1_data = A2_pos[:, 2]

                elif param == "vel":
                    A1_data = A2_vel

                elif param == "v":
                    A1_data = root_sum_sq(A2_vel.T)

                elif param == "v_x":
                    A1_data = A2_vel[:, 0]

                elif param == "v_y":
                    A1_data = A2_vel[:, 1]

                elif param == "v_z":
                    A1_data = A2_vel[:, 2]

                elif param == "rho":
                    A1_data = A1_rho

                elif param == "u":
                    A1_data = A1_u

                elif param == "phi":
                    A1_data = ut.load_particle_array(f, "phi")

                    # Account for an external potential
                    if self.impact is not None and self.impact.m_xp is not None:
                        A1_r_xp = root_sum_sq((A2_pos - self.impact.A1_pos_xp).T)
                        A1_data -= G * self.impact.m_xp / A1_r_xp

                elif param == "T":
                    A1_data = woma.A1_T_u_rho(A1_u, A1_rho, A1_mat_id % ut.id_body)

                elif param == "s":
                    A1_data = woma.A1_s_u_rho(A1_u, A1_rho, A1_mat_id % ut.id_body)

                elif param == "o":
                    A1_data = oset

                elif param == "a":
                    A1_data = oset.A1_a

                elif param == "e":
                    A1_data = oset.A1_e

                elif param == "i":
                    A1_data = oset.A1_i

                elif param == "q":
                    A1_data = oset.A1_q

                elif param == "Q":
                    A1_data = oset.A1_Q

                elif param == "a_eq":
                    A1_data = oset.A1_a * (1 - oset.A1_e**2) * np.cos(oset.A1_i) ** 2

                elif param == "Q_v_inf":
                    # Initialise all apoapses
                    A1_data = oset.A1_Q

                    # Overwrite with speeds at infinity for unbound particles
                    A1_sel_unb = np.where(oset.A1_e >= 1)[0]
                    A1_data[A1_sel_unb] = oset.A1_v_inf[A1_sel_unb]

                elif param == "fof_id":
                    # Load the fof-output snapshot file
                    Fp_fof_snap = Fp_snap[:-5] + "_fof_%s_0000.hdf5" % self.link_len

                    print('Loading "%s" ' % Fp_fof_snap[-52:])
                    with h5py.File(Fp_fof_snap, "r") as f_fof:
                        A1_fof_id = ut.load_particle_array(f_fof, param)

                        # Sort by the fof-file particle IDs
                        A1_id_fof = ut.load_particle_array(f_fof, "id")
                        A1_sort_fof = np.argsort(A1_id_fof)

                    A1_data = A1_fof_id[A1_sort_fof]

                else:
                    A1_data = ut.load_particle_array(f, param)[A1_sort]

                A1_return.append(A1_data)

        if len(A1_param) == 1:
            return A1_return[0]
        else:
            return A1_return

    def load_fof_data(self, time, A1_param, link_len=None, fof_id_max=None):
        """Load (or compute) one or more arrays of fof group data.

        Compute either the centre-of-mass, mean, or summed parameter value.

        Parameters
        ----------
        time : int
            The time of the snapshot to load (s).
            0       The initial conditions.
            -1      The final snapshot.

        A1_param : [str]
            A list of fof group properties to load, or one.

        link_len : str (opt.)
            The fof linking length (file_to_SI.l).

        fof_id_max : int (opt.)
            Only include fof groups up to and including this ID.

        Returns
        -------
        A1_return : [[?]]
            A list of the arrays of fof group properties, or one.
        """
        self.init()
        time = int(time)

        A1_return = []
        A1_param = ensure_list(A1_param)

        # Defaults
        if link_len is None:
            link_len = self.link_len
        if fof_id_max is None:
            fof_id_max = self.po.fof_id_max

        # Sum or mean parameters (otherwise centre of mass)
        A1_param_sum = ["m"]
        A1_param_mean = ["rho", "u", "T", "s"]

        # Load particle fof IDs and/or fof-group orbits if needed
        if any(
            param in A1_param_sum + A1_param_mean + ["fof_id"] for param in A1_param
        ):
            A1_fof_id = self.load_snapshot_data(time, "fof_id", link_len)
        if any(param not in A1_param_sum + A1_param_mean for param in A1_param):
            oset_fof = self.load_or_compute_oset_fof(time, link_len, fof_id_max)

        # ========
        # Load or compute each requested array
        # ========
        for param in A1_param:
            # Sum
            if param in A1_param_sum:
                A1_m_particle = self.load_snapshot_data(time, param)
                A1_data = [
                    sum(A1_m_particle[A1_fof_id == fof_id])
                    for fof_id in range(fof_id_max + 1)
                ]

            elif param == "f_c":
                A1_m_particle, A1_mat_id_particle = self.load_snapshot_data(
                    time, ["m", "mat_id"]
                )
                A1_sel_c = A1_mat_id_particle % ut.id_body == self.mat_id_c
                A1_data = [
                    sum(A1_m_particle[A1_sel_c & A1_fof_id == fof_id])
                    for fof_id in range(fof_id_max + 1)
                ]

            # Mean
            elif param in A1_param_mean:
                A1_p_particle = self.load_snapshot_data(time, param)
                A1_data = [
                    mean(A1_p_particle[A1_fof_id == fof_id])
                    for fof_id in range(fof_id_max + 1)
                ]

            # Centre of mass
            elif param == "pos":
                A1_data = oset_fof.A2_pos

            elif param == "r":
                A1_data = root_sum_sq(oset_fof.A2_pos.T)

            elif param == "x":
                A1_data = oset_fof.A2_pos[:, 0]

            elif param == "y":
                A1_data = oset_fof.A2_pos[:, 1]

            elif param == "z":
                A1_data = oset_fof.A2_pos[:, 2]

            elif param == "vel":
                A1_data = oset_fof.A2_vel

            elif param == "v":
                A1_data = root_sum_sq(oset_fof.A2_vel.T)

            elif param == "v_x":
                A1_data = oset_fof.A2_vel[:, 0]

            elif param == "v_y":
                A1_data = oset_fof.A2_vel[:, 1]

            elif param == "v_z":
                A1_data = oset_fof.A2_vel[:, 2]

            elif param == "a":
                A1_data = oset_fof.A1_a

            elif param == "e":
                A1_data = oset_fof.A1_e

            elif param == "i":
                A1_data = oset_fof.A1_i

            elif param == "q":
                A1_data = oset_fof.A1_q

            elif param == "Q":
                A1_data = oset_fof.A1_Q

            elif param == "a_eq":
                A1_data = (
                    oset_fof.A1_a
                    * (1 - oset_fof.A1_e**2)
                    * np.cos(oset_fof.A1_i) ** 2
                )

            elif param == "Q_v_inf":
                # Initialise all apoapses
                A1_data = oset_fof.A1_Q

                # Overwrite with speeds at infinity for unbound groups
                A1_sel_unb = np.where(oset_fof.A1_e >= 1)[0]
                A1_data[A1_sel_unb] = oset_fof.A1_v_inf[A1_sel_unb]

            # Misc
            elif param == "fof_id":
                A1_data = A1_fof_id

            else:
                raise Exception("Not set to load fof param ", param)

            A1_return.append(np.array(A1_data))

        if len(A1_param) == 1:
            return A1_return[0]
        else:
            return A1_return

    def param_value(self, param):
        """Extract a parameter from inputs or results, in plotting units."""
        self.init()

        # ========
        # Results
        # ========
        # Load pre-recorded values if they exist, or calculate and save them
        if param in ["m_capt", "m_c_capt"]:
            if "A2000c30" in self.name:
                (
                    m_bnd_in_Hill,
                    m_bnd_in_05Hill,
                    m_bnd_out_Hill,
                    m_unb,
                    m_bnd_in_Hill_fof_min,
                    num_fof_min,
                    num_fof_min_capt,
                    m_c_bnd_in_Hill,
                    m_c_bnd_in_05Hill,
                    m_c_bnd_out_Hill,
                    m_c_unb,
                ) = ut_ph.accum_data_phodei_m_capt(self)
            else:
                (
                    m_bnd_in_Hill,
                    m_bnd_in_05Hill,
                    m_bnd_out_Hill,
                    m_unb,
                    m_bnd_in_Hill_fof_min,
                    num_fof_min,
                    num_fof_min_capt,
                ) = ut_ph.accum_data_phodei_m_capt(self)

            if param == "m_capt":
                return m_bnd_in_Hill / self.impact.ic_2.M
            elif param == "m_c_capt":
                return m_c_bnd_in_Hill / self.impact.ic_2.M

        # ========
        # Inputs
        # ========
        # Custom options and/or units
        if param == "N":
            return self.num_picle
        elif param in ["L_1", "L_2"]:
            if param == "L_1":
                ic = self.impact.ic_1
            elif param == "L_2":
                ic = self.impact.ic_2

                # Manually distinguish differentiated asteroids
                if self.category == "phodei" and "A2000c30" in self.name:
                    return -1

            if ic.ip.period is None:
                return 0
            else:
                if hasattr(ic.ip, "L_max"):
                    return round_to_nearest(ic.ip.L / ic.ip.L_max, 0.125)
                else:
                    return ic.ip.L
        elif param in ["L_1_ax", "L_2_ax"]:
            if param == "L_1_ax":
                A1_preset = self.impact.A1_preset_1
            elif param == "L_2_ax":
                A1_preset = self.impact.A1_preset_2

            if "spin=x" in A1_preset:
                return 1
            elif "spin=y" in A1_preset:
                return 2
            elif "spin=z" in A1_preset:
                return 3
            elif "spin=-x" in A1_preset:
                return -1
            elif "spin=-y" in A1_preset:
                return -2
            elif "spin=-z" in A1_preset:
                return -3
            else:
                return 0

        # Evaluate directly, standard/preset units
        if param in dir(self):
            value = float(getattr(self, param))
        elif param in dir(self.impact):
            value = float(getattr(self.impact, param))

        unit = self.po.Di_param_unit_label[param].unit

        return value / unit

    # ========
    # Misc utilities
    # ========
    def select_snapshot_times(self, time=None, time_end=None, time_step=None):
        """Select a list of available snapshot times.

        Parameters
        ----------
        time : int (opt.)
            The time of the first (or only) snapshot to select (s).
            None        All snapshots.
            0           The initial conditions.
            -1          The final snapshot.
            "list=*"    Instead, load a list of times, where * is the file path.

        time_end : int (opt.)
            The time of the last snapshot to select, if time is not None.
            None        Only select one snapshot.
            -1          The final snapshot.

        time_step : int (opt.)
            The time steps between selected snapshots, if time_end is not None.
            None        All snapshots between time and time_end.

        Returns
        -------
        A1_time : [int]
            Snapshot times (s).
        """
        self.init()
        # List provided instead
        if time[:5] == "list=":
            Fp_snap_list = check_end(time[5:], ".txt")
            A1_time = np.loadtxt(Fp_snap_list).astype(float)

            # Snapshots
            A1_sel = np.in1d(self.A1_time, A1_time)
            A1_time = self.A1_time[A1_sel]

            if len(A1_time) == 0:
                raise Exception("No snapshots selected from %s" % Fp_snap_list)

            return A1_time

        time = check_int_or_def(time, None)
        time_end = check_int_or_def(time_end, None)
        time_step = check_int_or_def(time_step, None)

        # All
        if time is None:
            return self.A1_time

        # Convert if needed
        time = int(time)
        if time == 0:
            time = self.A1_time[0]
        elif time == -1:
            time = self.A1_time[-1]

        # Just one
        if time_end is None:
            if time in self.A1_time:
                return [time]
            else:
                raise Exception("No snapshot selected for %s at %d" % (self.name, time))

        # Convert if needed
        time_end = int(time_end)
        if time_end == -1:
            time_end = self.A1_time[-1]

        # All in range
        if time_step is None:
            A1_sel = np.where((self.A1_time > time) & (self.A1_time < time_end))[0]
            if len(A1_sel) > 0:
                return self.A1_time[A1_sel]
            else:
                raise Exception(
                    "No snapshots selected for %s between %d %d"
                    % (self.name, time, time_end)
                )

        # Convert if needed
        time_step = int(time_step)

        # Requested times
        A1_time = np.arange(time, time_end, time_step)

        # Snapshots
        A1_sel = np.in1d(np.round(self.A1_time, 5), np.round(A1_time, 5))
        if len(A1_sel) > 0:
            return self.A1_time[A1_sel]
        else:
            raise Exception(
                "No snapshots selected for %s in range %d, %d, %d"
                % (self.name, time, time_end, time_step)
            )

    def Fp_snap_from_time(self, time):
        """Return the corresponding snapshot file path for the given time.

        Parameters
        ----------
        time : int
            The time of the snapshot to load (s).
            0       The initial conditions.
            -1      The final snapshot.
        """
        time = int(time)

        if time == 0:
            return self.Fp_init_cond
        elif time == -1:
            return self.A1_Fp_snap[-1]

        idx_snap = idx_closest(self.A1_time, time)
        if abs(self.A1_time[idx_snap] - time) > 0.05:
            print(
                "\nWarning: closest snapshot loaded at %.1f not %.1f"
                % (self.A1_time[idx_snap], time)
            )

        return self.A1_Fp_snap[idx_snap]

    def init_orbit(self):
        """Compute the initial orbit of body 2 around body 1."""
        self.init()
        self.impact.prep_scenario()

        return Orbit(
            self.impact.A1_pos_2 - self.impact.A1_pos_1,
            self.impact.A1_vel_2 - self.impact.A1_vel_1,
            self.impact.M_1,
            self.impact.M_2,
        )

    def init_orbits(self):
        """Compute the initial orbits of both bodies around an external potential.

        Also return the centre-of-mass orbit.
        """
        self.init()

        # Get positions and velocities from stored rebound-collision info
        A1_pos_c_1 = self.Di_misc["A1_pos_c_1"]
        A1_vel_c_1 = self.Di_misc["A1_vel_c_1"]
        m_c_1 = self.Di_misc["m_c_1"]
        A1_pos_c_2 = self.Di_misc["A1_pos_c_2"]
        A1_vel_c_2 = self.Di_misc["A1_vel_c_2"]
        m_c_2 = self.Di_misc["m_c_2"]

        o_0_1 = Orbit(A1_pos_c_1, A1_vel_c_1, self.impact.m_xp, m_c_1)
        o_0_2 = Orbit(A1_pos_c_2, A1_vel_c_2, self.impact.m_xp, m_c_2)

        # Centre of mass
        o_0_c = Orbit(
            (m_c_1 * A1_pos_c_1 + m_c_2 * A1_pos_c_2) / (m_c_1 + m_c_2),
            (m_c_1 * A1_vel_c_1 + m_c_2 * A1_vel_c_2) / (m_c_1 + m_c_2),
            self.impact.m_xp,
            m_c_1 + m_c_2,
        )

        return o_0_1, o_0_2, o_0_c

    # ========
    # Plotting utilities
    # ========
    def tstamp(self, time):
        """Format the time stamp as a string for plotting annotations etc."""
        self.init()

        # Use the time to contact or periapsis for t=0
        if self.impact is not None:
            self.impact.prep_scenario()
            if self.impact.t_c is not None:
                time -= self.impact.t_c
            if self.impact.t_q is not None:
                time -= self.impact.t_q

        tstamp = round(time * s_to_hour, 1)

        # Avoid "-0"
        if abs(tstamp) <= 0.05:
            tstamp = 0.0

        return r"$%.1f$ h" % tstamp

    def param_value_str(self, param, value=None):
        """Convert a parameter value to a string for annotations etc.

        Parameters
        ----------
        param : str
            The parameter to convert, e.g. a simulation property.

        value : ? (opt.)
            If provided, then instead of extracting the value automatically from
            param_value, convert the provided value instead.
        """
        self.init()

        # Custom options
        if param == "N":
            N_log = np.log10(self.num_picle if value is None else value)

            # Round to nearest 0.25
            N_log = 0.25 * np.round(N_log / 0.25)

            return "$10^{%.2g}$" % N_log
        if param in ["L_1_ax", "L_2_ax"]:
            if value is None:
                value = self.param_value(param)

            # Decode spin axis
            if value == 0:
                return r"$0$"
            elif value == 1:
                return r"$x$"
            elif value == 2:
                return r"$y$"
            elif value == 3:
                return r"$z$"
            elif value == -1:
                return r"$-x$"
            elif value == -2:
                return r"$-y$"
            elif value == -3:
                return r"$-z$"

        # Convert directly to string
        else:
            return "%s" % (self.param_value(param) if value is None else value)

    def param_label_unit(self, param):
        """A label with units for latex plotting etc.

        Usually just a wrapper for Di_param_unit_label, but with custom options.
        """
        self.init()

        # Custom options
        if param in ["L_1", "L_2"]:
            if param == "L_1":
                ic = self.impact.ic_1
            elif param == "L_2":
                ic = self.impact.ic_2

            if hasattr(ic.ip, "L_max"):
                if self.name == "Ma_xp_A2000_n65_r11_v00":
                    return r"$L_z$ ($L_{\rm max}$)"
                else:
                    return r"Spin AM ($L_{\rm max}$)"
            else:
                return r"Spin AM (%s)" % Di_param_unit_label["L"].unit_label
        elif param in ["L_1_ax", "L_2_ax"]:
            return "Spin axis"
        elif param == "m_capt":
            return "Mass Fraction"

        # Standard/preset options
        else:
            return self.po.Di_param_unit_label[param].label_unit

    def add_text(self, time, po=None, ax=None):
        """Add text to a figure.

        Parameters
        ----------
        time : int
            The time of the snapshot being plotted.

        po : PlotOptions (opt.)
            Use instead of self.po, if provided.

        ax : Axes (opt.)
            The plot axes.
        """
        if po is None:
            po = self.po
        if ax is None:
            ax = plt.gca()

        # Axis limits
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        # Plot each bit of text
        for i_text, (loc, text) in enumerate(
            [
                ["upper left", po.text_ul],
                ["upper center", po.text_uc],
                ["upper right", po.text_ur],
                ["lower left", po.text_ll],
                ["lower center", po.text_lc],
                ["lower right", po.text_lr],
            ]
        ):
            if text is None:
                continue

            # Scale bar
            if text[:6] == "scale=":
                # Scale length
                if text[6:] == "auto":
                    scale_max = (x_max - x_min) / 4

                    # Order of magnitude
                    log10 = np.floor(np.log10(scale_max))
                    magn = 10**log10

                    # 1, 2, or 5 * 10^log
                    if scale_max / magn < 2:
                        scale = magn
                    elif scale_max / magn < 5:
                        scale = 2 * magn
                    else:
                        scale = 5 * magn
                else:
                    scale = float(text[6:])

                scalebar = AnchoredSizeBar(
                    transform=ax.transData,
                    size=scale,
                    label=r"$%g$ (%s)"
                    % (scale, po.Di_param_unit_label["r"].unit_label),
                    loc=loc,
                    color=po.text_colour,
                    sep=po.fontsize_text / 2,
                    pad=0.2,
                    frameon=False,
                    label_top="upper" in loc,
                    fontproperties=mpl.font_manager.FontProperties(
                        size=po.fontsize_text
                    ),
                )
                ax.add_artist(scalebar)

            # Text
            else:
                # Automatic labels
                if text == "tstamp":
                    text = self.tstamp(time)
                elif text[:5] == "eval=":
                    text = "%s" % eval(text[5:])

                add_text_loc(
                    text,
                    loc,
                    c=po.text_colour,
                    outl_c=po.rgba_bg,
                    fontsize=po.fontsize_text,
                    xy_ratio=(x_max - x_min) / (y_max - y_min),
                    ax=ax,
                )

    def plot_planet(self, A1_ax_lim=None, po=None, ax=None):
        """Plot external-potential planet (if within the axis limits).

        Also plot the Roche limit and/or Hill sphere if applicable.

        Parameters
        ----------
        A1_ax_lim : [float] (opt.)
            The axis limits: [x_0, x_1, y_0, y_1] (r_unit).

        po : PlotOptions (opt.)
            An object with settings for plotting functions.

        ax : Axes (opt.)
            The axes to plot on.

        Required po.Di_misc_colour
        --------------------------
        "planet" : float
            The colour for the planet.

        "Roche" : float (opt.)
            If provided, then the colour for the Roche limit.

        "Hill" : float (opt.)
            If provided, then the colour for the Hill sphere.
        """
        if po is None:
            po = self.po
        if ax is None:
            ax = plt.gca()

        # Plot if it fits within the axes
        R_planet = self.impact.ic_1.R_s / po.r_unit
        if A1_ax_lim is None or (
            A1_ax_lim[0] < R_planet
            and A1_ax_lim[1] > -R_planet
            and A1_ax_lim[2] < R_planet
            and A1_ax_lim[3] > -R_planet
        ):
            circ = plt.Circle(
                (0, 0),
                R_planet,
                facecolor=po.Di_misc_colour["planet"],
                edgecolor="none",
                alpha=0.7,
                zorder=-999,
            )
            ax.add_patch(circ)

            # Also plot the Roche limit
            if "Roche" in po.Di_misc_colour.keys():
                circ = plt.Circle(
                    (0, 0),
                    po.R_Roche / po.r_unit,
                    facecolor="none",
                    edgecolor=po.Di_misc_colour["Roche"],
                    lw=1,
                    ls=ls_dot,
                    alpha=0.4,
                    zorder=-999,
                )
                ax.add_patch(circ)

        # Colour outside of the Hill sphere
        if "Hill" in po.Di_misc_colour.keys():
            R_Hill = po.R_Hill / po.r_unit
            # Plot if visible
            if A1_ax_lim is None or any(A1_ax_lim / np.sqrt(2) > R_Hill):
                # Colour the whole axes, then plot a background-colour disk
                ax.set_facecolor(po.Di_misc_colour["Hill"])
                circ = plt.Circle(
                    (0, 0),
                    R_Hill,
                    facecolor=po.rgba_bg,
                    edgecolor="none",
                    zorder=-9999,
                )
                ax.add_patch(circ)

    # ========
    # Shortcut properties
    # ========
    @property
    def m_picle(self):
        """Inherit the mean particle mass from either the impact or ic object."""
        self.init()

        if self.impact is not None:
            return self.impact.m_picle
        else:
            return self.ic.m_picle

    @property
    def mat_id_c(self):
        """The core material ID."""
        self.init()

        if self.impact is not None:
            return self.impact.ic_2.A1_mat_id_layer[0]
        else:
            return self.ic.A1_mat_id_layer[0]

    @property
    def A1_pos_p(self):
        """The position of the primary for computing orbits, from an external potential."""
        return self.impact.A1_pos_xp

    @property
    def A1_vel_p(self):
        """The position of the primary for computing orbits, from an external potential."""
        return np.array([0, 0, 0])

    @property
    def M_p(self):
        """The position of the primary for computing orbits, from an external potential."""
        return self.impact.m_xp


class SimSet:
    """A set of Simulation objects, with things like plotting-option overrides.

    Most of the list attributes, apart from A1_sim, can be given as a single
    [value] to be automatically duplicated into [value] * len(A1_sim).

    Parameters
    ----------
    name : str
        The simulation set name.

    copy : str
        The name of a different SimSet object from which to copy all attributes
        that are not provided as other arguments here.

    A1_sim : [Simulation] (or [RebSim])
        The list of Simulation objects.

    A1_time : [int]
        A list of times to select a single snapshot from each simulation.
        0           The initial conditions.
        -1          The final snapshot.

    A1_label : [str]
        A list of labels for plotting. Defaults to ["sim *"] order numbers. Or
        set ["auto=*"] for some automatic label options.

    A1_colour : [str]
        A list of colours for plotting. Defaults to A1_c in jkeger.py. Set to
        [cmap] to automatically use evenly spaced colours from that colour map.

    A1_linestyle : [str]
        A list of line styles for plotting. Defaults to all solid lines.

    A1_marker : [str]
        A list of markers for plotting. Defaults to all circles.

    A1_param_c : [str]
        A list of particle properties to colour them by.

    A1_zorder : [str]
        A list of custom zorders.

    A2_bonus : [[str]]
        A list of bonus arguments to override ut.A1_bonus.

    po : PlotOptions
        An object with settings for plotting functions. Defaults to inherit from
        A1_sim[0] if not provided.

    A1_text_ul, _uc, _ur, _ll, _lc, _lr : [str]
        A list of text labels to override PlotOptions.text_*.

    legend_title : str
        A title for the legend.

    set_label : str
        A label for the set as a whole, to distinguish between other sets.

    set_label_title : str
        A title for a legend across other sets, along with `set_label`.

    Di_param_Di_p_colour, Di_param_Di_p_linestyle : {}
        Dictionaries of, for subset parameter names, a dictionary of subset
        parameter values and corresponding colours or linestyles.

    category : str
        The base simulation type. Defaults to A1_sim[0]'s category.

    do_share_y : bool
        Whether to share axes for tile plots.

    wspace : float
        Horizontal spacing between tile plots.

    no_min_x_tl, no_min_y_tl : bool
        If True then remove the lowest x, y tick label, e.g. to help fit with
        another row of panels underneath.

    Di_misc : {}
        A dictionary of miscellaneous parameters.

    Attributes
    ----------
    num_sim : int
        The number of simulations in the set.
    """

    # Dictionary of all SimSet objects
    _Di_sim_set = {}

    def __init__(
        self,
        name,
        copy=None,
        A1_sim=None,
        A1_time=None,
        A1_label=None,
        A1_colour=None,
        A1_linestyle=None,
        A1_marker=None,
        A1_param_c=None,
        A1_zorder=None,
        A2_bonus=None,
        po=None,
        A1_text_ul=None,
        A1_text_uc=None,
        A1_text_ur=None,
        A1_text_ll=None,
        A1_text_lc=None,
        A1_text_lr=None,
        legend_title=None,
        category=None,
        set_label=None,
        set_label_title=None,
        Di_param_Di_p_colour=None,
        Di_param_Di_p_linestyle=None,
        do_share_y=None,
        wspace=None,
        no_min_x_tl=None,
        no_min_y_tl=None,
        Di_misc=None,
    ):
        self.name = name
        self.copy = copy
        self.A1_sim = np.array(A1_sim)
        self.A1_time = A1_time
        self.A1_label = A1_label
        self.A1_colour = A1_colour
        self.A1_linestyle = A1_linestyle
        self.A1_marker = A1_marker
        self.A1_param_c = A1_param_c
        self.A1_zorder = A1_zorder
        self.A2_bonus = A2_bonus
        self.po = po
        self.A1_text_ul = A1_text_ul
        self.A1_text_uc = A1_text_uc
        self.A1_text_ur = A1_text_ur
        self.A1_text_ll = A1_text_ll
        self.A1_text_lc = A1_text_lc
        self.A1_text_lr = A1_text_lr
        self.legend_title = legend_title
        self.category = category
        self.set_label = set_label
        self.set_label_title = set_label_title
        self.Di_param_Di_p_colour = Di_param_Di_p_colour
        self.Di_param_Di_p_linestyle = Di_param_Di_p_linestyle
        self.do_share_y = do_share_y
        self.wspace = wspace
        self.no_min_x_tl = no_min_x_tl
        self.no_min_y_tl = no_min_y_tl
        self.Di_misc = Di_misc

        self.init_done = False

        # Add to the dictionary
        if self.name is not None:
            SimSet._Di_sim_set[self.name] = self

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        for sim in self.A1_sim:
            sim.init()
        self.num_sim = len(self.A1_sim)

        # Copy
        if self.copy is not None:
            copy = SimSet._Di_sim_set[self.copy]
            copy.init()
            copy_object(self, copy, A1_set_eq=["A1_sim"])

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        # Defaults
        if self.A1_label is None:
            self.A1_label = ["sim %d" for i in range(self.num_sim)]
        if self.A1_colour is None:
            self.A1_colour = [A1_c[i % len(A1_c)] for i in range(self.num_sim)]
        if self.A1_linestyle is None:
            self.A1_linestyle = ["-"] * self.num_sim
        if self.A1_marker is None:
            self.A1_marker = ["o"] * self.num_sim
        if self.po is None:
            self.po = self.A1_sim[0].po
        if self.category is None:
            self.category = self.A1_sim[0].category
        if self.do_share_y is None:
            self.do_share_y = True
        if self.no_min_x_tl is None:
            self.no_min_x_tl = False
        if self.no_min_y_tl is None:
            self.no_min_y_tl = False

        # Duplicate single inputs for all simulations
        if self.A1_time is not None and len(self.A1_time) != self.num_sim:
            self.A1_time = [self.A1_time[0]] * self.num_sim
        if len(self.A1_label) != self.num_sim:
            if self.A1_label[0][:5] == "auto=":
                auto_param = self.A1_label[0][5:]

                # Labels
                self.A1_label = [sim.param_value_str(auto_param) for sim in self.A1_sim]

                # Legend title
                if auto_param == "N":
                    self.legend_title = r"$N$"
                elif auto_param in sim.po.Di_param_unit_label:
                    self.legend_title = sim.po.Di_param_unit_label[auto_param].label
                else:
                    self.legend_title = None
            else:
                self.A1_label = [self.A1_label[0]] * self.num_sim
        if len(self.A1_colour) != self.num_sim:
            # Use colour map if provided
            try:
                self.A1_colour = self.A1_colour[0](np.linspace(0, 1, self.num_sim))
            except TypeError:
                self.A1_colour = [self.A1_colour[0]] * self.num_sim
        if len(self.A1_linestyle) != self.num_sim:
            self.A1_linestyle = [self.A1_linestyle[0]] * self.num_sim
        if len(self.A1_marker) != self.num_sim:
            self.A1_marker = [self.A1_marker[0]] * self.num_sim
        if self.A1_param_c is not None and len(self.A1_param_c) != self.num_sim:
            self.A1_param_c = [self.A1_param_c[0]] * self.num_sim
        if self.A2_bonus is not None and len(self.A2_bonus) != self.num_sim:
            self.A2_bonus = [self.A2_bonus[0]] * self.num_sim
        if self.A1_text_ul is not None and len(self.A1_text_ul) != self.num_sim:
            self.A1_text_ul = [self.A1_text_ul[0]] * self.num_sim
        if self.A1_text_uc is not None and len(self.A1_text_uc) != self.num_sim:
            self.A1_text_uc = [self.A1_text_uc[0]] * self.num_sim
        if self.A1_text_ur is not None and len(self.A1_text_ur) != self.num_sim:
            self.A1_text_ur = [self.A1_text_ur[0]] * self.num_sim
        if self.A1_text_ll is not None and len(self.A1_text_ll) != self.num_sim:
            self.A1_text_ll = [self.A1_text_ll[0]] * self.num_sim
        if self.A1_text_lc is not None and len(self.A1_text_lc) != self.num_sim:
            self.A1_text_lc = [self.A1_text_lc[0]] * self.num_sim
        if self.A1_text_lr is not None and len(self.A1_text_lr) != self.num_sim:
            self.A1_text_lr = [self.A1_text_lr[0]] * self.num_sim

        # Ensure array types
        self.A1_time = np.array(self.A1_time)
        self.A1_label = np.array(self.A1_label)
        self.A1_colour = np.array(self.A1_colour)
        self.A1_linestyle = np.array(self.A1_linestyle)
        self.A1_marker = np.array(self.A1_marker)
        self.A1_param_c = np.array(self.A1_param_c)
        self.A2_bonus = np.array(self.A2_bonus)
        self.A1_text_ul = np.array(self.A1_text_ul)
        self.A1_text_ur = np.array(self.A1_text_ur)
        self.A1_text_ll = np.array(self.A1_text_ll)
        self.A1_text_lr = np.array(self.A1_text_lr)

    def ensure_unique(self):
        """Ensure each simulation only appears once in the set."""
        A1_sim_unique = []
        A1_done = []
        for sim in self.A1_sim:
            if sim.name in A1_done:
                continue
            else:
                A1_done.append(sim.name)
                A1_sim_unique.append(sim)

        self.A1_sim = A1_sim_unique


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR classes_sim.py  ====\n")
