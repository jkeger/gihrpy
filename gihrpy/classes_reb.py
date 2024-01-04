"""
More classes for organising GIHR projects, specifically for rebound simulations.

See classes.py.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.utilities_reb as ur
from gihrpy.classes import PlotOptions


class RebPicle:
    """Input information for setting up a REBOUND particle.

    Parameters
    ----------
    m : float
        The mass (kg).

    R : float (opt.)
        The radius (m).

    A1_pos, A1_vel : [float] (opt.)
        The position and velocity (m, m/s).

    a, e, i, nu, Omega, pomega, omega : float (opt.)
         The semi-major axis, eccentricity, inclination, true anomaly, longitude
         of ascending node, longitude of periapsis, and argument of periapsis (SI).

     name : str (opt.)
         A name for reference.
    """

    def __init__(
        self,
        m,
        R=None,
        A1_pos=None,
        A1_vel=None,
        a=None,
        e=None,
        i=None,
        nu=None,
        Omega=None,
        pomega=None,
        omega=None,
        name=None,
    ):
        self.m = m
        self.R = R
        self.A1_pos = None if A1_pos is None else np.array(A1_pos).astype(float)
        self.A1_vel = None if A1_vel is None else np.array(A1_vel).astype(float)
        self.a = a
        self.e = e
        self.i = 0 if i is None else i
        self.nu = nu
        self.Omega = 0 if Omega is None else Omega
        self.pomega = 0 if pomega is None else pomega
        self.omega = 0 if omega is None else omega
        self.name = name


class RebSim:
    """Information for a REBOUND simulation.

    Parameters
    ----------
    name : str
        The simulation name.

    copy : str
        The name of a different RebSim object from which to copy all attributes
        that are not provided as other arguments here.

    po : PlotOptions
        An object with settings for plotting functions. Defaults to inherit from
        ic if not provided.

    dir_proj : str
        The file path to the base directory for the project, expected to contain
        a sub directory for this simulation.

    category : str
        The base type of simulation, e.g. the project name.

    t_start, t_out_step, t_end : float
        The start, output-step, and end times (yr).

    t_wall_restart : float
        The wall-clock time between writing restart files (h).

    t_2(_3) : float
        The time to change to the second (third) output step (s).

    t_out_step_2(_3) : float
        The second (third) time between outputs, from t_2 (t_3) onwards (s).

    Fp_init_cond : str
        The file path for the initial conditions. Defaults to
        <dir_proj>/<name>/init_cond.txt.

    Fp_config : str
        The file path for the config parameters. Defaults to
        <dir_proj>/<name>/config.txt.

    Fp_sim : str
        The file path for this simulation's source code. Defaults to
        <dir_proj>/<name>/simulation.c.

    A1_picle : [RebPicle]
        Input information for the particles.

    id_cent : int
        The particle ID for the central body for estimating orbits etc.

    id_orb_start : int
        The particle ID for the first particle orbiting the central body.

    sim_in : Simulation
        A Simulation object to use for the initial conditions.

    time_in : int
        The snapshot time (s) for the initial conditions from sim_in.

    num_fof : int
        The number of friends of friends groups to add as orbiting particles.

    m_fof_min : float
        The minimum mass (kg) of FoF groups to add, overrides num_fof.

    bnd_orb_Q_max : float
        If provided, then only add bound orbiting particles with Q below this.

    test_masses : str
        Options to set orbiting particles to be zero-mass test particles:

        "none"  (Default) No test particles.
        "all"   All orbiting particles as test masses.
        "m<*"   Orbiting particles below the given mass as test masses.

    sim_in_rot_ax, sim_in_angle : str, float
        The axis and angle (degrees) by which to rotate the input simulation.

    sim_in_rot_ax_2, sim_in_angle_2 : str, float
        The axis and angle (degrees) for a second rotation.

    n_split : int
        Split test particles between this many runs as crude parallelisation.

    id_rel_prim : int
        The ID of the designated primary particle to remove orbiting particles
        that get too close or too far from it.

    rel_distance_max : double
        The maximum distance from the designated primary particle beyond which
        particles will be removed.

    oblate_id : int
        The ID of the oblate particle.

    oblate_J2 : double
        The J2 moment of the oblate particle.

    oblate_obliquity : double
        The obliquity (rad) of the oblate particle.

    collision_mode : int
        The type of inter-particle collisions to use: -1 = no collisions;
        0 = print collision info only; 1 = hard-sphere collisions; 2 = merge.

    coeff_restitution : double
        The coefficient of restitution for hard-sphere collisions.

    t_collision_delay : double
        Delay collisions to only be allowed after this time has elapsed (s),
        e.g. to allow particles to separate from close initial conditions.

    radius_scale : double
        For collision detection, multiply up the particle radii by this value.

    id_frame_shift : int
        The ID of a particle to which the reference frame should be shifted
        every heartbeat. Can yield a big speed-up for e.g. a planet with
        orbiting particles that are together all orbiting a star.

    Di_po_edit : {}
        A dictionary of parameter names and values to edit the PlotOptions
        object.

    Di_misc : {}
        A dictionary of miscellaneous parameters.

    Attributes
    ----------
    num_picle : int
        The total number of particles.

    num_test : int
        The number of test particles.

    A1_time : [int]
        A list of all snapshot times (yr).

    A1_snap_id : [str]
        The list of snapshot IDs as strings.

    Fp_output : [str]
        The file path for the output data.

    Fp_collisions : [str]
        The file path for the collisions data.

    Fp_log : [str]
        The file path for the log.

    A1_m : [float]
        The masses of the particles (kg). May not be complete at first.

    A1_R : [float]
        The radii of the particles (m). May not be complete at first.

    A2_pos : [[float]]
        The positions of the particles (m). May not be complete at first.

    A2_vel : [[float]]
        The velocities of the particles (m/s). May not be complete at first.

    Fp_accum_data : str
        The file path for accumulated derived-output etc data.
    """

    # Dictionary of all RebSim objects
    _Di_reb_sim = {}

    def __init__(
        self,
        name,
        copy=None,
        po=None,
        dir_proj=None,
        category=None,
        t_start=None,
        t_out_step=None,
        t_end=None,
        t_wall_restart=None,
        t_2=None,
        t_out_step_2=None,
        t_3=None,
        t_out_step_3=None,
        Fp_init_cond=None,
        Fp_config=None,
        Fp_sim=None,
        A1_picle=None,
        id_cent=None,
        id_orb_start=None,
        sim_in=None,
        time_in=None,
        num_fof=None,
        m_fof_min=None,
        bnd_orb_Q_max=None,
        test_masses=None,
        sim_in_rot_ax=None,
        sim_in_angle=None,
        sim_in_rot_ax_2=None,
        sim_in_angle_2=None,
        n_split=None,
        id_rel_prim=None,
        rel_distance_max=None,
        oblate_id=None,
        oblate_J2=None,
        oblate_obliquity=None,
        collision_mode=None,
        coeff_restitution=None,
        t_collision_delay=None,
        radius_scale=None,
        id_frame_shift=None,
        Di_po_edit=None,
        Di_misc=None,
    ):
        self.name = name
        self.copy = copy
        self.po = po
        self.dir_proj = dir_proj
        self.category = category
        self.t_start = t_start
        self.t_out_step = t_out_step
        self.t_end = t_end
        self.t_wall_restart = t_wall_restart
        self.t_2 = t_2
        self.t_out_step_2 = t_out_step_2
        self.t_3 = t_3
        self.t_out_step_3 = t_out_step_3
        self.Fp_init_cond = Fp_init_cond
        self.Fp_config = Fp_config
        self.Fp_sim = Fp_sim
        self.A1_picle = A1_picle
        self.id_cent = id_cent
        self.id_orb_start = id_orb_start
        self.sim_in = sim_in
        self.time_in = time_in
        self.num_fof = num_fof
        self.m_fof_min = m_fof_min
        self.bnd_orb_Q_max = bnd_orb_Q_max
        self.test_masses = test_masses
        self.sim_in_rot_ax = sim_in_rot_ax
        self.sim_in_angle = sim_in_angle
        self.sim_in_rot_ax_2 = sim_in_rot_ax_2
        self.sim_in_angle_2 = sim_in_angle_2
        self.n_split = n_split
        self.id_rel_prim = id_rel_prim
        self.rel_distance_max = rel_distance_max
        self.oblate_id = oblate_id
        self.oblate_J2 = oblate_J2
        self.oblate_obliquity = oblate_obliquity
        self.collision_mode = collision_mode
        self.coeff_restitution = coeff_restitution
        self.t_collision_delay = t_collision_delay
        self.radius_scale = radius_scale
        self.id_frame_shift = id_frame_shift
        self.Di_po_edit = Di_po_edit
        self.Di_misc = Di_misc

        self.init_done = False

        # Add to the dictionary
        if self.name is not None:
            RebSim._Di_reb_sim[self.name] = self

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        # Copy
        if self.copy is not None:
            copy = RebSim._Di_reb_sim[self.copy]
            copy.init()
            copy_object(
                self,
                copy,
                A1_exclude=["Fp_init_cond", "Fp_config", "Fp_sim"],
                A1_set_eq=["po"],
            )

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        if self.sim_in is not None:
            self.sim_in.init()

        # Defaults
        if self.po is None:
            self.po = PlotOptions()
        if self.Fp_init_cond is None:
            self.Fp_init_cond = "%s/%s/init_cond.txt" % (self.dir_proj, self.name)
        if self.Fp_config is None:
            self.Fp_config = "%s/%s/config.txt" % (self.dir_proj, self.name)
        if self.Fp_sim is None:
            self.Fp_sim = "%s/%s/simulation.c" % (self.dir_proj, self.name)
        if self.t_start is None:
            self.t_start = 0
        if self.id_cent is None:
            self.id_cent = 0
        if self.id_orb_start is None:
            self.id_orb_start = 0
        if self.num_fof is None:
            self.num_fof = 0
        if self.test_masses is None:
            self.test_masses = "none"

        # Derived attributes
        self.Fp_stem = "%s/%s" % (self.dir_proj, self.name)
        self.A1_time = np.round(
            np.arange(self.t_start, self.t_end + 1e-6, self.t_out_step), 5
        )
        self.A1_snap_id = np.array(["%.2f" % time for time in self.A1_time])
        self.Fp_output = "%s/output.txt" % self.Fp_stem
        self.Fp_collisions = "%s/collisions.txt" % self.Fp_stem
        self.Fp_log = "%s/log.txt" % self.Fp_stem
        self.A1_m = np.array([p.m for p in self.A1_picle])
        self.A1_R = np.array([p.R for p in self.A1_picle])
        self.A2_pos = np.zeros((len(self.A1_picle), 3))
        self.A2_vel = np.zeros((len(self.A1_picle), 3))
        self.Fp_accum_data = "%s_accum_data.hdf5" % self.Fp_stem

        # Overrides etc
        if self.Di_po_edit is not None:
            # Make a clean copy to edit
            self.po = deepcopy(self.po)

            for var in self.Di_po_edit.keys():
                # self.po.var = self.Di_po_edit["var"]
                exec("self.po." + var + ' = self.Di_po_edit["' + var + '"]')

    def formatted_label(self):
        """Return a formatted label for e.g. figure titles."""
        if self.sim_in_rot_ax == "z":
            phi = self.sim_in_angle
        elif self.sim_in_rot_ax_2 == "z":
            phi = self.sim_in_angle_2
        if self.sim_in_rot_ax == "x":
            i = -self.sim_in_angle
        else:
            i = 0
        label = r"$\phi = %d^\circ$, $\; i = %d^\circ$" % (phi, i)

        if self.oblate_obliquity is not None:
            o = -self.oblate_obliquity
            label += r", $\; o_{\small \mars{}} = %d^\circ$" % o
        if "_nu" in self.name:
            label += r", $\; e_{\small \mars{}} = %.3g$" % self.A1_picle[self.id_cent].e
            label += r", $\; \nu_{\small \mars{}} = %d^\circ$" % (
                self.A1_picle[self.id_cent].nu * rad_to_deg
            )
        if "_noSu" in self.name:
            label += r" (no Sun)"

        return label

    # ========
    # Data loading etc
    # ========
    def reb_load_or_compute_oset(self, num_orb=None, num_time=None, do_recalc=False):
        """Compute and save the particle orbits, or load if previously saved.

        Index order: oset.A2_o[idx_time, idx_orb].

        Parameters
        ----------
        num_orb : int (opt.)
            The number of orbiting particles required, default all. Recalculate
            if a smaller number were previously saved.

        num_time : float (opt.)
            The number of time steps required, default all. Recalculate if a
            smaller number were previously saved.

        do_recalc : bool (opt.)
            Force recalculation even if a previously saved file exists.

        Returns
        -------
        oset : OrbitSet
            The set of particles orbits.
        """
        self.init()

        Fp_oset = "%s_oset.pkl" % (self.Fp_stem)

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
            # Load output data
            (
                A1_time,
                A1_m,
                A1_R,
                A3_pos,
                A3_vel,
                A1_idx_last,
                A1_fate,
            ) = self.reb_load_data()

            if num_orb is None:
                num_orb = len(A1_m) - self.id_orb_start

            # Extract particles and times to compute
            A3_pos = A3_pos[:num_time, self.id_orb_start : self.id_orb_start + num_orb]
            A3_vel = A3_vel[:num_time, self.id_orb_start : self.id_orb_start + num_orb]
            A2_M_p = np.full((num_time, num_orb), A1_m[self.id_cent])
            A2_m = np.full(
                (num_time, num_orb),
                A1_m[self.id_orb_start : self.id_orb_start + num_orb],
            )

            print("Computing %s orbits... " % (self.name[-48:]), flush=True, end="")
            oset = OrbitSet(A3_pos=A3_pos, A3_vel=A3_vel, A2_M_p=A2_M_p, A2_m=A2_m)
            print("Done")

            print('\rWriting "%s" ' % Fp_oset[-58:], flush=True, end="")
            with open(Fp_oset, "wb") as f_oset:
                pickle.dump(oset, f_oset)
            print("Done      ")
        else:
            print('Loading "%s" ' % Fp_oset[-58:], flush=True, end="")
            with open(Fp_oset, "rb") as f_oset:
                oset = pickle.load(f_oset)
            print("Done")

            # Force recalculation if the required particles and times are not in the saved set
            if (num_orb is not None and oset.num_orb_2 < num_orb) or (
                num_time is not None and oset.num_orb_1 < num_time
            ):
                print(
                    "Recompute... (num_orb %d, num_time %d)"
                    % (oset.num_orb_2, oset.num_orb_1)
                )
                self.reb_load_or_compute_oset(
                    num_orb=num_orb, num_time=num_time, do_recalc=True
                )

        return oset

    # Function wrappers
    def reb_load_data(self, A1_time_sel=None, do_reconvert=False):
        A1_time, A1_m, A1_R, A3_pos, A3_vel, A1_idx_last, A1_fate = ur.reb_load_data(
            self, A1_time_sel=A1_time_sel, do_reconvert=do_reconvert
        )

        return (
            A1_time,
            A1_m,
            A1_R,
            A3_pos,
            A3_vel,
            A1_idx_last,
            A1_fate,
        )

    def reb_load_snapshot(self, time):
        return ur.reb_load_snapshot(self, time)

    def reb_load_collisions(self, do_reconvert=False, do_no_repeats=True):
        return ur.reb_load_collisions(
            self, do_reconvert=do_reconvert, do_no_repeats=do_no_repeats
        )


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR classes_reb.py  ====\n")
