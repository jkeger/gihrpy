"""
General utilities for GIHR REBOUND projects.
"""

from jkeger import *

# Internal units (for G = 1)
reb_m_unit = 1.988475e30  # kg (M_Sun)
reb_l_unit = 1.4959787e11  # m (AU)
reb_t_unit = 5.02254803e6  # s (yr/2pi)

# Flags etc
reb_collision_mode_none = -1
reb_collision_mode_print = 0
reb_collision_mode_hard_sphere = 1
reb_collision_mode_merge = 2
reb_fate_survived = -1
reb_fate_rm_close = 1
reb_fate_rm_far = 2
reb_fate_rm_merge = 3


# ========
# Initial conditions
# ========
def reb_write_config(
    Fp_save,
    num_picle,
    num_test,
    file_init_cond,
    file_out_stem,
    file_restart,
    t_start,
    t_end,
    t_out_step,
    t_wall_restart=None,
    t_2=None,
    t_out_step_2=None,
    t_3=None,
    t_out_step_3=None,
    id_rel_prim=None,
    rel_distance_max=None,
    oblate_id=None,
    oblate_J2=None,
    oblate_obliquity=None,
    collision_mode=None,
    coeff_restitution=None,
    t_collision_delay=None,
    id_frame_shift=None,
):
    """Make a rebound config file.

    File contents (-1 for default values or to skip an optional parameter)
    -------------
    # HEADER (8 lines)
    num_picle  num_test                                     (-  -)
    file_init_cond  file_out_stem  file_restart             (-  -  -)
    t_start  t_end  t_out_step  t_wall_restart              (yr  yr  yr  h)
    t_2  t_out_step_2  t_3  t_out_step_3                    (yr  yr  yr  yr)
    id_rel_prim  rel_distance_max                           (-  m)
    oblate_id  oblate_J2  oblate_obliquity                  (-  -  rad)
    collision_mode  coeff_restitution  t_collision_delay    (-  -  yr)
    id_frame_shift                                          (-)                                             (-)

    Parameters
    ----------
    Fp_save : str
        The file path to write to.

    See RebSim for the rest.
    """
    print('Saving to "%s"... ' % Fp_save[-52:], end="")
    with open(Fp_save, "w") as f:
        # Header
        f.write("# num_picle  num_test                                     (-  -) \n")
        f.write(
            "# file_init_cond  file_out_stem  file_restart             (-  -  -) \n"
        )
        f.write(
            "# t_start  t_end  t_out_step  t_wall_restart              (yr  yr  yr  h) \n"
        )
        f.write(
            "# t_2  t_out_step_2  t_3  t_out_step_3                    (yr  yr  yr  yr) \n"
        )
        f.write("# id_rel_prim  rel_distance_max                           (-  m) \n")
        f.write(
            "# oblate_id  oblate_J2  oblate_obliquity                  (-  -  rad) \n"
        )
        f.write(
            "# collision_mode  coeff_restitution  t_collision_delay    (-  -  yr) \n"
        )
        f.write("# id_frame_shift                                          (-) \n")

        # Default/skip values
        if t_wall_restart is None:
            t_wall_restart = 10.0  # h
        if t_2 is None:
            t_2 = -1.0
            t_out_step_2 = -1.0
        if t_3 is None:
            t_3 = -1.0
            t_out_step_3 = -1.0
        if id_rel_prim is None:
            id_rel_prim = -1
            rel_distance_max = -1.0
        if oblate_id is None:
            oblate_id = -1
            oblate_J2 = -1.0
            oblate_obliquity = -1.0
        if collision_mode is None:
            collision_mode = -1
        if coeff_restitution is None:
            coeff_restitution = -1.0
        if t_collision_delay is None:
            t_collision_delay = -1.0
        if id_frame_shift is None:
            id_frame_shift = -1

        # Data
        f.write("%d %d \n" % (num_picle, num_test))
        f.write("%s %s %s \n" % (file_init_cond, file_out_stem, file_restart))
        f.write("%.8e %.8e %.8e %.3f \n" % (t_start, t_end, t_out_step, t_wall_restart))
        f.write("%.8e %.8e %.8e %.8e \n" % (t_2, t_out_step_2, t_3, t_out_step_3))
        f.write("%d %.8e \n" % (id_rel_prim, rel_distance_max))
        f.write(
            "%d %.8e %.8e \n" % (oblate_id, oblate_J2, oblate_obliquity * deg_to_rad)
        )
        f.write(
            "%d %.8e %.8e \n" % (collision_mode, coeff_restitution, t_collision_delay)
        )
        f.write("%d \n" % (id_frame_shift))

    print("Done")


def reb_write_init_cond(Fp_save, num_picle, A1_m, A1_R, A2_pos, A2_vel):
    """Make a rebound initial conditions file. See RebSim.

    File contents
    -------------
    # HEADER (2 lines)
    num_picle                                       (-)
    m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0    (kg  m  m  m  m  m/s  m/s  m/s)
    m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                 ...                  ...
    m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n

    Parameters
    ----------
    Fp_save : str
        The file path to write to.

    See RebSim for the rest.
    """
    print('Saving to "%s"... ' % Fp_save[-52:], end="")
    # Save as a combined array
    np.savetxt(
        Fp_save,
        np.concatenate(([A1_m], [A1_R], A2_pos.T, A2_vel.T)).T,
        fmt=["%.8e", "%.8e", "%.12e", "%.12e", "%.12e", "%.12e", "%.12e", "%.12e"],
        header=(
            "# num_picle                        (-) \n"
            "# m  R  x  y  z  v_x  v_y  v_z     (kg  m  m  m  m  m/s  m/s  m/s) \n"
            "%d" % num_picle
        ),
        comments="",
    )
    print("Done")


# ========
# Loading
# ========
def reb_load_data(rsim, A1_time_sel=None, do_reconvert=False):
    """Load particle data from a rebound output file.

    Returned positions and velocities are relative to the central mass, if any.

    File contents
    -------------
    # t  m  R  x  y  z  v_x  v_y  v_z    (yr  kg  m  m  m  m  m/s  m/s  m/s)
    t_0  m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0
    t_0  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                    ...                    ...
    t_0  m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n
    t_1  m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0
    t_1  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                    ...                    ...
    t_1  m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n
    ...                    ...                    ...

    Parameters
    ----------
    rsim : RebSim
        The rebound simulation object.

    A1_time_sel : [float] (opt.)
        Selected snapshot times to extract (yr). Defaults to all.

    do_reconvert : bool (opt.)
        Set True to overwrite any converted data and load from the original.
        Automatically set True if the output file was more recently modified.

    Returns
    -------
    A1_time : [float]
        The output times of each snapshot (yr).

    A1_m : [float]
        The particle masses (kg).

    A1_R : [float]
        The particle radii (m).

    A3_pos : [[[float]]]
        The particle positions, at each time (m).

    A3_vel : [[[float]]]
        The particle velocities, at each time (m s^-1).

    A1_idx_last : [float]
        The index of the last snapshot each particle was in before removal.

    A1_fate : [int]
        The fate of each particle:
            -1      Survived
            1       Too close to primary
            2       Too far from primary
            3       Merge collision
    """
    print("Loading %s data... " % rsim.name[-42:], end="", flush=True)

    # Input and converted data file paths (for split 0)
    if rsim.n_split is None or rsim.n_split < 2:
        Fp_output = rsim.Fp_output
        rsim.n_split = 0
    else:
        Fp_output = rsim.Fp_output[:-4] + "_0.txt"
    Fp_output_conv = rsim.Fp_output[:-4] + ".npz"

    # Placeholder simulations with no orbiting particles
    if rsim.name in [
        "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
    ]:
        return (
            np.array([0, 5e3]),
            np.array([Su.M, Ma.M]),
            np.array([Su.R, Ma.R]),
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [Ma.a, 0.0, 0.0],
                ]
            ),
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.0, Ma.v_orb, 0.0],
                ]
            ),
            np.ones(2).astype(int) * -1,
            np.ones(2).astype(int) * reb_fate_survived,
        )

    # Check if the data file was modified more recently than converted data
    if (
        not os.path.isfile(Fp_output_conv)
        or os.path.getmtime(Fp_output) > os.path.getmtime(Fp_output_conv) + 3600
    ):
        do_reconvert = True

    # Load converted numpy data if it exists
    if not do_reconvert and os.path.isfile(Fp_output_conv):
        with np.load(Fp_output_conv) as A2_data:
            A1_time_data = A2_data["t"]
            A1_m_data = A2_data["m"]
            A1_R_data = A2_data["R"]
            A2_pos_data = A2_data["pos"]
            A2_vel_data = A2_data["vel"]

        # No need to load separate splits, should already be included
        do_load_splits = False
    else:
        do_reconvert = True

        # Load the output data
        # A2_data = np.loadtxt(Fp_output, skiprows=1)
        A2_data = np.genfromtxt(Fp_output, skip_header=1, invalid_raise=False)

        # Remove data with tiny time differences
        dt_min = 1e-3
        A1_t = np.unique(A2_data[:, 0])
        A1_dt = np.diff(A1_t)
        A1_t_sel = A1_t[np.where(A1_dt > dt_min)[0]]
        A1_sel = np.isin(A2_data[:, 0], A1_t_sel)
        A2_data = A2_data[A1_sel]

        # Parse the data
        A1_time_data = np.copy(A2_data[:, 0])
        A1_m_data = np.copy(A2_data[:, 1])
        A1_R_data = np.copy(A2_data[:, 2])
        A2_pos_data = np.copy(A2_data[:, 3:6])
        A2_vel_data = np.copy(A2_data[:, 6:])

        # Append any split-run data
        for i_split in range(1, rsim.n_split):
            # File path for this split
            Fp_output = rsim.Fp_output[:-4] + "_%d.txt" % i_split

            # Load the output data
            A2_data = np.loadtxt(Fp_output, skiprows=1)

            # Remove data with tiny time differences
            A1_t = np.unique(A2_data[:, 0])
            A1_dt = np.diff(A1_t)
            A1_t_sel = A1_t[np.where(A1_dt > dt_min)[0]]
            A1_sel = np.isin(A2_data[:, 0], A1_t_sel)
            A2_data = A2_data[A1_sel]

            # Parse the split data
            A1_time_data_i = np.copy(A2_data[:, 0])
            A1_m_data_i = np.copy(A2_data[:, 1])
            A1_R_data_i = np.copy(A2_data[:, 2])
            A2_pos_data_i = np.copy(A2_data[:, 3:6])
            A2_vel_data_i = np.copy(A2_data[:, 6:])

            # Overwrite the new time values to equal the closest existing match
            A1_time_0 = np.unique(A1_time_data)
            A1_time_i = np.unique(A1_time_data_i)
            for time_0 in A1_time_0:
                time_i = A1_time_i[idx_closest(A1_time_i, time_0)]
                A1_time_data_i[A1_time_data_i == time_i] = time_0

            # Only append matched times ## or later times?
            A1_sel_append = np.where(
                np.isin(
                    A1_time_data_i, A1_time_0
                )  # | (A1_time_data_i > np.max(A1_time_0))
            )

            # Append
            A1_time_data = np.concatenate((A1_time_data, A1_time_data_i[A1_sel_append]))
            A1_m_data = np.concatenate((A1_m_data, A1_m_data_i[A1_sel_append]))
            A1_R_data = np.concatenate((A1_R_data, A1_R_data_i[A1_sel_append]))
            A2_pos_data = np.concatenate((A2_pos_data, A2_pos_data_i[A1_sel_append]))
            A2_vel_data = np.concatenate((A2_vel_data, A2_vel_data_i[A1_sel_append]))

    # Save numpy data
    if do_reconvert:
        np.savez(
            Fp_output_conv,
            t=A1_time_data,
            m=A1_m_data,
            R=A1_R_data,
            pos=A2_pos_data,
            vel=A2_vel_data,
        )

    # Snapshot times
    A1_time = np.unique(A1_time_data)

    # Extract a subset of times
    if A1_time_sel is not None:
        # Take each closest time
        A1_time = [A1_time[idx_closest(A1_time, time_sel)] for time_sel in A1_time_sel]

    # Masses and radii from the initial conditions
    Fp_init_cond = check_end(rsim.Fp_init_cond, ".txt")
    A1_m = np.loadtxt(Fp_init_cond, skiprows=3, usecols=(0))
    A1_R = np.loadtxt(Fp_init_cond, skiprows=3, usecols=(1))
    num_picle = len(A1_R)

    # Initialise the arrays
    A3_pos = np.full((len(A1_time), num_picle, 3), np.nan)
    A3_vel = np.full((len(A1_time), num_picle, 3), np.nan)
    A2_pos = np.full((num_picle, 3), np.nan)
    A2_vel = np.full((num_picle, 3), np.nan)

    # Set each snapshot time's data
    for i_time, time in enumerate(A1_time):
        # Select this time
        A1_sel_time = np.where(A1_time_data == time)[0]

        # Set each particle's data by matching radii, nan for removed particles
        A2_pos.fill(np.nan)
        A2_vel.fill(np.nan)
        for i_picle_snap, R in enumerate(A1_R_data[A1_sel_time]):
            i_picle = idx_closest(A1_R, R)
            A2_pos[i_picle] = A2_pos_data[A1_sel_time][i_picle_snap]
            A2_vel[i_picle] = A2_vel_data[A1_sel_time][i_picle_snap]

        # Shift to the central mass frame
        if rsim.id_cent is not None:
            A2_pos -= A2_pos[rsim.id_cent]
            A2_vel -= A2_vel[rsim.id_cent]

        # Set data
        A3_pos[i_time] = A2_pos
        A3_vel[i_time] = A2_vel

    # Load log file info for particle removal comments
    ## allow for multiple log files?
    try:
        #   Header (~15 lines after the start of the file):
        # Steps      Particles   Time (yr) (frac)   Timestep (yr) CPU time (s) (tot, min)
        #   Example standard and removal lines:
        # 3900348    N=50        t=24.000007 ( 0.5%)  dt=0.000025  wall=3.106 (27.7)
        # # Remove particle 27, r=1.39777918e-07: d=2.27635e-05 < 2.26566e-05 (t=31.1)
        # # Remove particle 79, r=9.38356021e-08: merge with 9, r=4.11094951e-07 (t=0.64)
        A2_log = np.genfromtxt(
            rsim.Fp_log,
            skip_header=17,
            skip_footer=9,
            comments=None,
            dtype=str,
            usecols=(4, 6),  # ["r=1234e5:", "<" or ">" or "with"] (on removal lines)
            invalid_raise=False,
        )

        A1_R_rm = np.copy(A2_log[:, 0])
        A1_sel_rm = np.char.startswith(A1_R_rm, "r=")
        A1_R_rm = np.char.strip(A1_R_rm[A1_sel_rm], chars="r=:").astype(float) * au_to_m
        A1_ineq_rm = np.copy(A2_log[:, 1])[A1_sel_rm]
    except IndexError:
        A1_R_rm = []
        A1_ineq_rm = []

    # Removed particles
    A1_idx_last = np.ones(num_picle).astype(int) * -1
    A1_fate = np.ones(num_picle).astype(int) * reb_fate_survived

    # Removal type
    for i_picle in range(num_picle):
        A1_x = A3_pos[:, i_picle, 0]
        if np.isnan(A1_x[-1]):
            # Find the first non-nan element in the reversed-order x values
            idx_last = len(A1_x) - (np.argmax(~np.isnan(A1_x)[::-1]) + 1)
            A1_idx_last[i_picle] = idx_last

            # Removal type from log file, by matching radii
            idx_rm = idx_closest(A1_R_rm, A1_R[i_picle])
            if A1_ineq_rm[idx_rm] == "<":
                A1_fate[i_picle] = reb_fate_rm_close
            elif A1_ineq_rm[idx_rm] == ">":
                A1_fate[i_picle] = reb_fate_rm_far
            elif A1_ineq_rm[idx_rm] == "with":
                A1_fate[i_picle] = reb_fate_rm_merge
            else:
                raise Exception(
                    "Failed to parse particle removal log:\n  ", A1_ineq_rm[idx_rm]
                )

    print("Done")

    return A1_time, A1_m, A1_R, A3_pos, A3_vel, A1_idx_last, A1_fate


def reb_load_snapshot(rsim, time):
    """Load particle data at a specific time.

    Parameters
    ----------
    rsim : RebSim
        The rebound simulation object.

    time : float
        The time to extract (yr).

    Returns
    -------
    A1_m : [float]
        The particle masses (kg).

    A1_R : [float]
        The particle radii (m).

    A2_pos : [[float]]
        The particle positions (m), relative to the central mass, if any.

    A2_vel : [[float]]
        The particle velocities (m s^-1), relative to the central mass, if any.
    """
    # Load all data
    A1_time, A1_m, A1_R, A3_pos, A3_vel, A1_idx_last, oset = reb_load_data(rsim)

    # Select this time (or the closest to it)
    if abs(A1_time[idx_time] - time) / min(A1_time[idx_time], time) > 0.5:
        raise Exception(
            "Closest snapshot at %.2f not %.2f yr" % (A1_time[idx_time], time)
        )
    elif abs(A1_time[idx_time] - time) > 0.05:
        print(
            "\nWarning: closest snapshot loaded at %.2f not %.2f yr"
            % (A1_time[idx_time], time)
        )

    return A1_m, A1_R, A3_pos[idx_time], A3_vel[idx_time]


def reb_load_collisions(rsim, do_reconvert=False, do_no_repeats=True):
    """Load collisions data from a rebound file.

    Assume no split runs for now.

    Returned positions and velocities are relative to the central mass, if any.

    File contents
    -------------
    # HEADER (1 line)
    t  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1  <same for _2>  r_prim  v_rel  (yr  SI...)

    Parameters
    ----------
    rsim : RebSim
        The rebound simulation object.

    do_reconvert : bool (opt.)
        Set True to overwrite any converted data and load from the original.
        Automatically set True if the output file was more recently modified.

    do_no_repeats : bool (opt.)
        If True, then discard collisions between the same pair of particles.

    Returns
    -------
    A1_time_c : [float]
        The time of each collision (yr).

    A1_r_prim_c : [float]
        The radial distance from the primary (m).

    A1_v_rel_c : [float]
        The relative velocity (m s^-1).

    A1_id_1_c, A1_id_2_c : [float]
        The IDs of the colliding particles, where id_1 has the greater mass.

    A2_pos_1_c, A2_vel_1_c, A2_pos_2_c, A2_vel_2_c : [float]
        The position (m) and velocity (m s^-1) of the colliding particles.
    """
    # No collisions
    if not os.path.isfile(rsim.Fp_collisions):
        return None, None, None, None, None, None, None, None, None

    print("Loading %s collision data... " % rsim.name[-42:], end="", flush=True)

    Fp_collisions_conv = rsim.Fp_collisions[:-4] + ".npz"

    # Check if the data file was modified more recently than converted data or a code change
    mtime_code_update = datetime.datetime.timestamp(
        datetime.datetime(2023, 9, 10, 1, 50)
    )
    if (
        not os.path.isfile(Fp_collisions_conv)
        or os.path.getmtime(rsim.Fp_collisions)
        > os.path.getmtime(Fp_collisions_conv) + 3600
        or os.path.getmtime(Fp_collisions_conv) < mtime_code_update
    ):
        do_reconvert = True

    # Load converted numpy data if it exists
    if not do_reconvert and os.path.isfile(Fp_collisions_conv):
        with np.load(Fp_collisions_conv) as A2_data:
            A1_time_c = A2_data["time"]
            A1_r_prim_c = A2_data["r_prim"]
            A1_v_rel_c = A2_data["v_rel"]
            A1_id_1_c = A2_data["id_1"]
            A1_id_2_c = A2_data["id_2"]
            A2_pos_1_c = A2_data["A2_pos_1"]
            A2_vel_1_c = A2_data["A2_vel_1"]
            A2_pos_2_c = A2_data["A2_pos_2"]
            A2_vel_2_c = A2_data["A2_vel_2"]
    else:
        # Load the data
        A2_data = np.loadtxt(rsim.Fp_collisions, skiprows=1)

        # No or invalid collisions data
        try:
            A2_data[:, 0]
        except IndexError:
            print("\n# Failed to parse collision data for", rsim.name)
            return None, None, None, None, None, None, None, None, None

        # Parse the data
        # t  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1  <same for *_2>  r_prim  v_rel
        # (yr) (SI)...
        A1_time_c = np.copy(A2_data[:, 0])
        offset_2 = 8
        A1_R_1_c = np.copy(A2_data[:, 2])
        A2_pos_1_c = np.copy(A2_data[:, 3:6])
        A2_vel_1_c = np.copy(A2_data[:, 6:9])
        A1_R_2_c = np.copy(A2_data[:, offset_2 + 2])
        A2_pos_2_c = np.copy(A2_data[:, offset_2 + 3 : offset_2 + 6])
        A2_vel_2_c = np.copy(A2_data[:, offset_2 + 6 : offset_2 + 9])
        A1_r_prim_c = np.copy(A2_data[:, 17])
        A1_v_rel_c = np.copy(A2_data[:, 18])
        num_c = len(A1_time_c)

        # Masses and radii from the initial conditions
        Fp_init_cond = check_end(rsim.Fp_init_cond, ".txt")
        A1_m = np.loadtxt(Fp_init_cond, skiprows=3, usecols=(0))
        A1_R = np.loadtxt(Fp_init_cond, skiprows=3, usecols=(1))

        # Set particle masses and IDs by matching the radii
        A1_R_unique = np.unique(np.append(A1_R_1_c, A1_R_2_c))
        Di_R_id = {R: idx_closest(A1_R, R) for R in A1_R_unique}
        A1_m_1_c = np.zeros(num_c)
        A1_m_2_c = np.zeros(num_c)
        A1_id_1_c = np.zeros(num_c, dtype=int)
        A1_id_2_c = np.zeros(num_c, dtype=int)
        for R, id in Di_R_id.items():
            A1_m_1_c[A1_R_1_c == R] = A1_m[id]
            A1_m_2_c[A1_R_2_c == R] = A1_m[id]
            A1_id_1_c[A1_R_1_c == R] = id
            A1_id_2_c[A1_R_2_c == R] = id
        # Set 1 particles to have the higher mass
        A1_sel_swap = np.where(A1_m_1_c < A1_m_2_c)[0]
        tmp = A1_id_1_c[A1_sel_swap]
        A1_id_1_c[A1_sel_swap] = A1_id_2_c[A1_sel_swap]
        A1_id_2_c[A1_sel_swap] = tmp

        # Save numpy data
        np.savez(
            Fp_collisions_conv,
            time=A1_time_c,
            r_prim=A1_r_prim_c,
            v_rel=A1_v_rel_c,
            id_1=A1_id_1_c,
            id_2=A1_id_2_c,
            A2_pos_1=A2_pos_1_c,
            A2_vel_1=A2_vel_1_c,
            A2_pos_2=A2_pos_2_c,
            A2_vel_2=A2_vel_2_c,
        )

    # Discard collisions between the same pair of particles
    if do_no_repeats:
        A1_pair = []
        A1_i_c_remove = []
        # Record duplicate pairs
        for i_c in range(len(A1_time_c)):
            pair = (A1_id_1_c[i_c], A1_id_2_c[i_c])
            if pair in A1_pair:
                A1_i_c_remove.append(i_c)
            else:
                A1_pair.append(pair)

            # Also remove unwantedly-recorded collisions with the primary
            if A1_id_1_c[i_c] == rsim.id_rel_prim:
                A1_i_c_remove.append(i_c)

        # Remove repeats
        A1_time_c = np.delete(A1_time_c, A1_i_c_remove)
        A1_r_prim_c = np.delete(A1_r_prim_c, A1_i_c_remove)
        A1_v_rel_c = np.delete(A1_v_rel_c, A1_i_c_remove)
        A1_id_1_c = np.delete(A1_id_1_c, A1_i_c_remove)
        A1_id_2_c = np.delete(A1_id_2_c, A1_i_c_remove)
        A2_pos_1_c = np.delete(A2_pos_1_c, A1_i_c_remove, axis=0)
        A2_vel_1_c = np.delete(A2_vel_1_c, A1_i_c_remove, axis=0)
        A2_pos_2_c = np.delete(A2_pos_2_c, A1_i_c_remove, axis=0)
        A2_vel_2_c = np.delete(A2_vel_2_c, A1_i_c_remove, axis=0)

    print("Done")

    return (
        A1_time_c,
        A1_r_prim_c,
        A1_v_rel_c,
        A1_id_1_c,
        A1_id_2_c,
        A2_pos_1_c,
        A2_vel_1_c,
        A2_pos_2_c,
        A2_vel_2_c,
    )


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR utilities_reb.py  ====\n")
