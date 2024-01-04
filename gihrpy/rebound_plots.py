"""
Functions for GIHR plots of REBOUND results, etc.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.utilities_reb as ur
from gihrpy.classes_sim import Simulation, SimSet
from gihrpy.classes_reb import RebSim
from gihrpy.objects import Di_sim, Di_sim_set, Di_reb_sim


def reb_plot_orbit_evol(
    rsim,
    A1_param_y=None,
    param_c=None,
    time_split=None,
    time_split_2=None,
    time_end=None,
    A1_idx_orb=None,
    dt=None,
    dt_1=None,
    dt_2=None,
):
    """Plot rebound particle estimated orbits through time over all snapshots.

    Plot multiple y parameters on separate x-aligned axes below the first.

    Parameters
    ----------
    rsim : RebSim or SimSet
        An object with all the information about the REBOUND simulation or a set
        of simulations.

    A1_param_y : [str] (opt.)
        The y parameter(s) to plot. Enter as a list including the brackets,
        commas, and no spaces, e.g. `[e,i,q]`.

    param_c : str (opt.)
        How to colour the particles (overridden by simulation sets), e.g. in
        order of: idx (default, affected by A1_idx_orb), m, e, Q, etc.

    time_split, time_split_2 : float (opt.)
        Split the x axis at this time (yr) to e.g. zoom in on early times, etc.

    time_end : float (opt.)
        An override end time (yr).

    A1_idx_orb : [int] (opt.)
        The indices of the orbiting particles to plot, default all. Enter as a
        list including the brackets, commas, and no spaces, e.g. `[0,4,7]`. Or
        enter e.g. `:5` to plot the first 5 particles (by default mass order).
        Or `custom` to extract from rsim.Di_misc.

    dt, dt_1, dt_2 : float (opt.)
        Sample the plotted times (in each split) at these intervals, e.g. for
        consistency across changes in the frequency of data outputs.

    Bonus options (set in ut.A1_bonus)
    -------------
    mark_rm         Indicate removal with a marker on the x axis.
    mark_c          Indicate collisions with a marker on the particle lines.
    """
    # Set the object if the name was provided instead
    if not isinstance(rsim, RebSim):
        try:
            rsim = Di_reb_sim[rsim]
            rsim.init()
            sim_set = SimSet(
                name="sim_set",
                A1_sim=[rsim],
            )
            sim_set.init()
        except KeyError:
            sim_set = Di_sim_set[rsim]
            sim_set.init()
            rsim = sim_set.A1_sim[0]
            rsim.init()
    num_sim = len(sim_set.A1_sim)

    A1_param_y = check_none_or_def(A1_param_y, "[e,i,q]")
    A1_param_y = array_from_string(A1_param_y, str)
    num_param_y = len(A1_param_y)
    param_c = check_none_or_def(param_c, "idx")
    time_split = check_float_or_def(time_split, None)
    time_split_2 = check_float_or_def(time_split_2, None)
    time_end = check_float_or_def(time_end, None)
    A1_idx_orb = check_none_or_def(A1_idx_orb, None)
    if A1_idx_orb == "custom":
        A1_idx_orb = rsim.Di_misc["A1_idx_orb"]
    elif A1_idx_orb is not None:
        if A1_idx_orb[0] == ":":
            A1_idx_orb = np.arange(int(A1_idx_orb[1:]))
        else:
            A1_idx_orb = array_from_string(A1_idx_orb, int)
        A1_idx_orb = np.sort(A1_idx_orb).astype(int)
    dt = check_float_or_def(dt, None)
    dt_1 = check_float_or_def(dt_1, None)
    dt_2 = check_float_or_def(dt_2, None)

    # Figure
    fig = plt.figure(
        figsize=(14 if num_param_y > 2 else 9, 9 + min(0, 3 * num_param_y - 3))
    )
    gs = mpl.gridspec.GridSpec(nrows=num_param_y, ncols=11)
    # Optionally split into columns
    if time_split_2 is not None:
        A1_ax_t0 = [plt.subplot(gs[i_ax, :2]) for i_ax in range(num_param_y)]
        A1_ax_t1 = [plt.subplot(gs[i_ax, 2:6]) for i_ax in range(num_param_y)]
        A1_ax_t2 = [plt.subplot(gs[i_ax, 6:]) for i_ax in range(num_param_y)]
    elif time_split is not None:
        A1_ax_t0 = [plt.subplot(gs[i_ax, :2]) for i_ax in range(num_param_y)]
        A1_ax_t1 = [plt.subplot(gs[i_ax, 2:]) for i_ax in range(num_param_y)]
        A1_ax_t2 = []
    else:
        A1_ax_t0 = [plt.subplot(gs[i_ax, :]) for i_ax in range(num_param_y)]
        A1_ax_t1 = []
        A1_ax_t2 = []
    gs.update(hspace=0.02, wspace=0)
    num_split = 1 + int(time_split is not None) + int(time_split_2 is not None)

    # Set a colour for each simulation in the set
    if num_sim > 1:
        A1_colour = cmap_rbow(np.linspace(0, 1, num_sim))
        param_c = "i_sim"

    # y limits for each axis, initialise to nan
    A2_y_lim = np.full((num_param_y, 2), np.nan)

    # ========
    # Plot each simulation
    # ========
    for i_rsim, rsim in enumerate(sim_set.A1_sim):
        po = rsim.po

        # Load output data
        A1_time, A1_m, A1_R, A3_pos, A3_vel, A1_idx_last, A1_fate = rsim.reb_load_data()
        A1_time_last = np.array(
            [A1_time[idx] if idx != -1 else -1 for idx in A1_idx_last]
        )

        # Store for axis limits
        time_start = np.round(A1_time[0], 0)

        # Select the particles to plot
        num_orb = len(A1_m) - rsim.id_orb_start
        if A1_idx_orb is not None:
            A1_idx_orb_plot = A1_idx_orb[A1_idx_orb < num_orb]
        else:
            A1_idx_orb_plot = np.arange(num_orb, dtype=int)
        num_orb_plot = len(A1_idx_orb_plot)

        # Override end time
        if time_end is None:
            time_end = A1_time[-1]
        A1_sel_time = A1_time <= time_end
        A1_time_last[A1_time_last >= time_end] = -1

        # Orbit data
        oset = rsim.reb_load_or_compute_oset(
            num_orb=np.amax(A1_idx_orb_plot) + 1, num_time=len(A1_time[A1_sel_time])
        )
        oset_sel = OrbitSet(A2_o=oset.A2_o[:, A1_idx_orb_plot])

        # Load collisions data
        if (
            any([param_y in A1_param_y for param_y in A1_param_y])
            or "mark_c" in ut.A1_bonus
        ):
            (
                A1_time_c,
                A1_r_prim_c,
                A1_v_rel_c,
                A1_id_1_c,
                A1_id_2_c,
                A2_pos_1_c,
                A2_vel_1_c,
                A2_pos_2_c,
                A2_vel_2_c,
            ) = rsim.reb_load_collisions()
            if A1_time_c is not None:
                A1_m_1_c = np.array([A1_m[id] for id in A1_id_1_c])
                A1_m_2_c = np.array([A1_m[id] for id in A1_id_2_c])
                A1_R_1_c = np.array([A1_R[id] for id in A1_id_1_c])
                A1_R_2_c = np.array([A1_R[id] for id in A1_id_2_c])

        # Colour each particle
        if num_sim == 1:
            # Avoid the too-light yellow end of plasma
            A1_colour = plt.get_cmap("plasma")(
                np.linspace(0, 0.97, len(A1_idx_orb_plot))
            )

        # ========
        # Plot
        # ========
        for i_split in range(num_split):
            # Select axes and times if splitting
            if time_split is None:
                A1_sel_time_split = np.ones(len(A1_time)).astype(bool)
                A1_ax = A1_ax_t0
            elif i_split == 0:
                A1_sel_time_split = (A1_time) <= (time_split + 1)
                A1_ax = A1_ax_t0
            elif i_split == 1:
                A1_sel_time_split = (A1_time >= time_split) & (
                    A1_time <= time_split_2 + 1
                )
                A1_ax = A1_ax_t1
            elif i_split == 2:
                A1_sel_time_split = A1_time >= time_split_2
                A1_ax = A1_ax_t2
            A1_sel_time_split = A1_sel_time_split & A1_sel_time

            # Extract a subset sampling of times
            dt_i = [dt, dt_1, dt_2][i_split]
            if dt_i is not None:
                # Round the times to their nearest dt
                A1_time_rounded = round_to_nearest(A1_time, dt_i)

                # Select only one time for each dt step
                A1_idx_sample = np.unique(A1_time_rounded, return_index=True)[1]
                A1_sel_sample = np.in1d(np.arange(len(A1_time)), A1_idx_sample)

                A1_sel_time_split = A1_sel_time_split & A1_sel_sample

            if len(A1_sel_time_split) == 0:
                continue
            time_split_start = A1_time[A1_sel_time_split][0]
            time_split_end = A1_time[A1_sel_time_split][-1]

            # ========
            # Plot each selected particle
            # ========
            for i_orb, idx_orb in enumerate(A1_idx_orb_plot):
                i_picle = rsim.id_orb_start + idx_orb
                # Note: oset_sel already has only selected particles, so is indexed with i_orb

                # Colour, linestyle, zorder
                if num_sim > 1:
                    # Colour by simulation
                    c = A1_colour[i_rsim]
                    # Linestyle by particle
                    ls = A1_ls[i_orb]
                    zorder = 2
                else:
                    if param_c == "idx":
                        # Colour by particle order
                        c = A1_colour[i_orb]
                        zorder = 2 + idx / num_orb_plot
                    else:
                        # Colour by particle initial-parameter order
                        A1_c_sort = np.argsort(oset_sel.A2_param(param_c))[0, :]
                        sel = np.where(A1_c_sort == i_orb)[0][0]
                        c = A1_colour[sel]
                        zorder = num_orb_plot - sel
                    ls = A1_ls[i_rsim]

                # Slightly thinner lines in later-time panels
                lw = [1.0, 0.9, 0.6][i_split]
                alpha = 0.93

                # ========
                # Plot each panel and parameter evolution
                # ========
                for i_ax, (ax, param_y) in enumerate(zip(A1_ax, A1_param_y)):
                    # Extract the parameter evolution
                    A1_p_y = oset_sel.A2_param(param_y)[:, i_orb]

                    # Convert to plotting units
                    A1_p_y /= po.Di_param_unit_label[param_y].unit

                    # Plot line
                    ax.plot(
                        A1_time[A1_sel_time_split],
                        A1_p_y[A1_sel_time_split],
                        c=c,
                        ls=ls,
                        lw=lw,
                        alpha=alpha,
                        zorder=zorder,
                    )

                    # Inspect late-time values for y limits
                    if i_split is None or i_split == 1:
                        A1_sel_time_late = A1_time >= 0.5 * time_end
                        y_min = np.nanmin(A1_p_y[A1_sel_time_late])
                        y_max = np.nanmax(A1_p_y[A1_sel_time_late])

                        # Store the accumulated min and max
                        if not np.isnan(y_min) and (
                            y_min < A2_y_lim[i_ax, 0] or np.isnan(A2_y_lim[i_ax, 0])
                        ):
                            A2_y_lim[i_ax, 0] = y_min
                        if not np.isnan(y_max) and (
                            y_max > A2_y_lim[i_ax, 1] or np.isnan(A2_y_lim[i_ax, 1])
                        ):
                            A2_y_lim[i_ax, 1] = y_max

                    # Label if in a set
                    if num_sim > 1 and i_split == 0 and i_ax == 0 and i_orb == 0:
                        ax.plot([], [], c=c, ls=ls, label=sim_set.A1_label[i_rsim])

                    # Mark collisions on the top panel
                    if "mark_c" in ut.A1_bonus:
                        A1_sel_c = np.where(
                            ((A1_id_1_c == i_picle) | (A1_id_2_c == i_picle))
                            & (A1_time_c > time_split_start)
                            & (A1_time_c < time_split_end)
                        )[0]

                        # Position on this particle's plotted line
                        A1_p_y_c = [
                            A1_p_y[idx_closest(A1_time, time_c)]
                            for time_c in A1_time_c[A1_sel_c]
                        ]

                        ax.scatter(
                            A1_time_c[A1_sel_c],
                            A1_p_y_c,
                            marker="*",
                            c=[c],
                            edgecolor="none",
                            s=6**2,
                            alpha=0.9,
                            zorder=2.5,
                        )

                # Mark particle removal times
                if "mark_rm" in ut.A1_bonus and A1_time_last[i_picle] != -1:
                    # Check the removal time is on this axis
                    x_rm = A1_time_last[i_picle] + rsim.t_out_step / 2
                    if x_rm > time_split_start and x_rm < time_split_end:
                        # Marker by type of removal
                        lw = 1
                        s = 5.5**2
                        if A1_fate[i_picle] == ur.reb_fate_rm_close:
                            marker = "x"
                            fc = [c]
                            ec = c
                        elif A1_fate[i_picle] == ur.reb_fate_rm_far:
                            marker = "o"
                            fc = "none"
                            ec = tweak_luminosity(c, 0.47)
                            lw = 1.7
                            s = 6**2
                        elif A1_fate[i_picle] == ur.reb_fate_rm_merge:
                            marker = "D"
                            fc = "none"
                            ec = c

                        # Position at bottom of y axis
                        y_rm = -2.1e-2
                        ax.scatter(
                            [x_rm],
                            [y_rm],
                            marker=marker,
                            lw=lw,
                            facecolor=fc,
                            edgecolor=ec,
                            s=s,
                            alpha=0.85,
                            transform=ax.get_xaxis_transform(),
                            zorder=2.5,
                            clip_on=False,
                        )

    # Legend
    if num_sim > 1:
        if time_split is not None:
            A1_ax = A1_ax_t0

        # Compact legend
        A1_ax[0].legend(
            prop={"size": 12},
            title=sim_set.legend_title,
            ncol=2,
            borderpad=0.3,
            handlelength=1,
            columnspacing=1,
        )

    # ========
    # Axes etc
    # ========
    # Labels
    if time_split is None:
        A1_ax_t0[-1].set_xlabel("Time (yr)")
    else:
        A1_ax_t1[-1].set_xlabel("Time (yr)")
    for i_ax, (ax, param_y) in enumerate(zip(A1_ax_t0, A1_param_y)):
        ax.set_ylabel(po.Di_param_unit_label[param_y].label_unit)

    # No tick labels on not-lowest rows
    for A1_ax in [A1_ax_t0, A1_ax_t1, A1_ax_t2]:
        for ax in A1_ax[:-1]:
            ax.set_xticklabels([])

    # x limits
    for ax in A1_ax_t0:
        ax.set_xlim(time_start, time_split)
    for ax in A1_ax_t1:
        ax.set_xlim(time_split, time_split_2)
    for ax in A1_ax_t2:
        ax.set_xlim(time_split_2, time_end)

    # y limits
    for i_ax, (ax, param_y) in enumerate(zip(A1_ax_t0, A1_param_y)):
        if param_y in po.Di_param_lim.keys():
            y_lim_0, y_lim_1 = po.Di_param_lim[param_y]
            ax.set_ylim(y_lim_0, y_lim_1)
        else:
            set_auto_limits(ax, [0, 0], A2_y_lim[i_ax], set_x=False, set_y=True)

    # Keep consistent y limits, and no ticks on not-leftmost columns
    for i_ax, ax in enumerate(A1_ax_t1):
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_ylim(A1_ax_t0[i_ax].get_ylim())
    for i_ax, ax in enumerate(A1_ax_t2):
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_ylim(A1_ax_t0[i_ax].get_ylim())

    # Misc extras
    if po.func_reb_orbit_evol_extra is not None:
        po.func_reb_orbit_evol_extra(
            rsim, A1_param_y, time_split, time_split_2, A1_ax_t0, A1_ax_t1, A1_ax_t2
        )

    A1_ax_all = A1_ax_t0 + A1_ax_t1 + A1_ax_t2
    for ax in A1_ax_all:
        auto_nice_plot(ax, fontsize=po.fontsize - 8, do_far_ticks=False)

    # Far axis ticks only on rightmost columns
    for i_ax, ax in enumerate(A1_ax_all[-num_param_y:]):
        ax_r = add_far_ticks(ax, do_x=False)[1]
        ax_r.set_yticks(A1_ax_t0[i_ax].get_yticks(minor=True), minor=True)
    tight_layout()

    # Save and close
    if num_sim == 1:
        Fp_save = "%s/%s/%s_evol" % (po.dir_save, rsim.name, rsim.name)
    else:
        Fp_save = "%s/%s_evol" % (po.dir_save, sim_set.name)
    for param_y in A1_param_y:
        Fp_save += "_%s" % param_y
    Fp_save += "_%s" % param_c
    if po.Fp_suffix is not None:
        Fp_save += "%s" % po.Fp_suffix

    # Save
    if ut.do_paper:
        savefig(Fp_save + ".pdf", dpi=400)
    else:
        savefig(Fp_save + ".png", dpi=po.dpi)
    print('\nSaved "%s" ' % Fp_save[-64:])

    plt.close()


def reb_plot_pos_trace(rsim, idx_orb, time=None, time_end=None):
    """Plot particle positions traced out through time.

    Parameters
    ----------
    rsim : RebSim
        An object with all the information about the REBOUND simulation.

    idx_orb : int
        The index of the orbiting particle to plot.

    time, time_end : float (opt.)
        The start and end times to plot between (yr).
        0       The initial simulation time.
        -1      The final simulation time.

    Bonus options (set in ut.A1_bonus)
    -------------
    xz_panel        Add an x--z projection panel above the main x--y one.
    zy_panel        Add a z--y projection panel right of the main x--y one.
    scale=*         Add a scale bar of this size (r_unit) instead of axis ticks,
                    set *=-1 for automatic.
    """
    # Set the object if the name was provided instead
    if not isinstance(rsim, RebSim):
        rsim = Di_reb_sim[rsim]
    rsim.init()
    po = rsim.po

    idx_orb = check_int(idx_orb)
    i_picle = rsim.id_orb_start + idx_orb
    time = check_float_or_def(time, 0)
    time_end = check_float_or_def(time_end, -1)

    # Load output data
    A1_time, A1_m, A1_R, A3_pos, A3_vel, A1_idx_last, A1_fate = rsim.reb_load_data()

    # Select time range
    if time == 0:
        time = A1_time[0]
    if time_end == -1:
        time_end = A1_time[-1]
    A1_time_sel = A1_time[np.where((A1_time >= time) & (A1_time <= time_end))[0]]
    if A1_idx_last[i_picle] != -1:
        i_time_last = np.where(A1_time_sel == A1_time[A1_idx_last[i_picle]])[0][0]
    else:
        i_time_last = len(A1_time_sel) - 2

    # Extract the selected particle's data
    oset = rsim.reb_load_or_compute_oset(num_orb=idx_orb + 1, num_time=i_time_last)
    oset = OrbitSet(A1_o=oset.A2_o[:, idx_orb])
    A2_pos = A3_pos[:, i_picle]
    A2_vel = A3_vel[:, i_picle]

    # Figure
    plt.figure()
    ax = plt.gca()
    ax.set_aspect("equal")

    # Add z-projection axes
    divider = make_axes_locatable(ax)
    divider.set_aspect(True)
    if "xz_panel" in ut.A1_bonus:
        if po.xz_size is not None:
            size = po.xz_size
        else:
            size = "50%"
        ax_xz = divider.append_axes("top", size=size, sharex=ax, pad=0)
        ax_xz.set_aspect("equal", adjustable="datalim")
    else:
        ax_xz = None
    if "zy_panel" in ut.A1_bonus:
        if po.zy_size is not None:
            size = po.zy_size
        else:
            size = "50%"
        ax_zy = divider.append_axes("top", size=size, sharex=ax, pad=0)
        ax_zy.set_aspect("equal", adjustable="datalim")
    else:
        ax_zy = None

    # Plot primary
    for ax_i in [ax, ax_xz, ax_zy]:
        if ax_i is None:
            continue

        circ = plt.Circle(
            (0, 0),
            A1_R[rsim.id_cent] / po.r_unit,
            facecolor="k",
            edgecolor="none",
            alpha=0.5,
        )
        ax_i.add_patch(circ)

    # ========
    # Interpolate particle positions
    # ========
    A2_pos_orig = np.copy(A2_pos)
    idx_insert = 0
    # Insert extra positions between each pair of times
    for i_time in range(i_time_last):
        # Index before which to insert extra points
        idx_insert += 1

        # Add additional positions using the predicted orbits either side
        o_0 = oset.A1_o[i_time]
        o_1 = oset.A1_o[i_time + 1]
        nu_0 = o_0.nu
        nu_1 = o_1.nu
        if nu_0 > nu_1:
            nu_1 += 2 * pi
        nu_step = nu_1 - nu_0

        # Add extra points to span between the angle gap
        # Finer steps around apoapsis where delta-nu means big position changes
        if o_0.e > 0.5 and abs(nu_0 * rad_to_deg - 180) < 10:
            nu_step_max = 1 * deg_to_rad
        elif o_0.e > 0.3 and abs(nu_0 * rad_to_deg - 180) < 30:
            nu_step_max = 5 * deg_to_rad
        elif o_0.e > 0.1 and abs(nu_0 * rad_to_deg - 180) < 120:
            nu_step_max = 10 * deg_to_rad
        else:
            nu_step_max = 30 * deg_to_rad
        n_intp = int(np.floor(nu_step / nu_step_max))

        # Interpolate the position at each nu from both orbits
        A1_nu_intp = np.linspace(nu_0, nu_1, n_intp + 2)
        A2_pos_intp = [
            (1 - i_nu / n_intp) * o_0.pos_vel_from_nu(nu, do_vel=False)
            + i_nu / n_intp * o_1.pos_vel_from_nu(nu, do_vel=False)
            for i_nu, nu in enumerate(A1_nu_intp)
        ][1:-1]

        # Insert into positions array
        A2_pos = np.insert(A2_pos, idx_insert, A2_pos_intp, axis=0)

        # Account for added elememts
        idx_insert += n_intp

    # ========
    # Plot
    # ========
    for ax_i, x_co, y_co in [[ax, 0, 1], [ax_xz, 0, 2], [ax_zy, 2, 1]]:
        if ax_i is None:
            continue

        plot_lines(
            A2_pos[: idx_insert + 1, x_co] / po.r_unit,
            A2_pos[: idx_insert + 1, y_co] / po.r_unit,
            A1_p_c=np.linspace(0, 1, idx_insert),
            cmap=mpl.colors.LinearSegmentedColormap.from_list(
                "", plt.get_cmap("viridis")(np.linspace(0, 0.82, 1024))
            ),
            lw=0.07 if ut.do_paper else 0.1,
            alpha=0.2 if ut.do_paper else 0.4,
            ax=ax_i,
        )

        # Mark if the particle was removed (within the plotted time)
        if A1_idx_last[i_picle] != -1 and i_time_last < len(A1_time_sel) - 1:
            ax_i.scatter(
                A2_pos_orig[i_time_last, x_co] / po.r_unit,
                A2_pos_orig[i_time_last, y_co] / po.r_unit,
                marker=".",
                c="r",
                s=3**2,
                edgecolor="none",
                alpha=0.7,
                zorder=99,
            )

    # Mark z=0 lines
    if ax_xz is not None:
        plt.setp(ax_xz.get_xticklabels(), visible=False)
        plt.setp(ax_xz.get_xticklines(), visible=False)
        ax_xz.axhline(0, c="k", lw=0.2, zorder=-9e9)
    if ax_zy is not None:
        plt.setp(ax_zy.get_yticklabels(), visible=False)
        plt.setp(ax_zy.get_yticklines(), visible=False)
        ax_zy.axvline(0, c="k", lw=0.2, zorder=-9e9)

    # Axis limits
    if "x" in po.Di_param_lim.keys():
        x_lim_0, x_lim_1 = po.Di_param_lim[param_y]
    else:
        x_min = np.nanmin(A2_pos[: idx_insert + 1, 0]) / po.r_unit
        x_max = np.nanmax(A2_pos[: idx_insert + 1, 0]) / po.r_unit
        x_lim_0 = x_min - (x_max - x_min) * 0.05
        x_lim_1 = x_max + (x_max - x_min) * 0.05
    ax.set_xlim(x_lim_0, x_lim_1)
    if "y" in po.Di_param_lim.keys():
        y_lim_0, y_lim_1 = po.Di_param_lim[param_y]
    else:
        y_min = np.nanmin(A2_pos[: idx_insert + 1, 1]) / po.r_unit
        y_max = np.nanmax(A2_pos[: idx_insert + 1, 1]) / po.r_unit
        y_lim_0 = y_min - (y_max - y_min) * 0.05
        y_lim_1 = y_max + (y_max - y_min) * 0.05
    ax.set_ylim(y_lim_0, y_lim_1)
    if ax_xz is not None or ax_zy is not None:
        if "z" in po.Di_param_lim.keys():
            z_lim_0, z_lim_1 = po.Di_param_lim[param_y]
        else:
            z_min = np.nanmin(A2_pos[: idx_insert + 1, 2]) / po.r_unit
            z_max = np.nanmax(A2_pos[: idx_insert + 1, 2]) / po.r_unit
            z_lim_0 = z_min - (z_max - z_min) * 0.05
            z_lim_1 = z_max + (z_max - z_min) * 0.05
        if ax_xz is not None:
            ax_xz.set_ylim(z_lim_0, z_lim_1)
        if ax_zy is not None:
            ax_zy.set_xlim(z_lim_0, z_lim_1)

    # Scale bar or axis ticks
    scale = ut.extract_bonus("scale")
    if scale is not None:
        # Auto scale
        if scale == "auto":
            scale_max = (ax.get_xlim()[1] - ax.get_xlim()[0]) / 4

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
            scale = float(scale)

        # Plot scale bar
        scalebar = AnchoredSizeBar(
            transform=ax.transData,
            size=scale,
            label=r"$%g$ %s" % (scale, po.r_unit_label),
            loc="lower left",
            color="k",
            sep=po.fontsize_text / 2,
            pad=0.2,
            frameon=False,
            label_top=False,
            fontproperties=mpl.font_manager.FontProperties(size=po.fontsize_text),
        )
        ax.add_artist(scalebar)

        # No axis ticks
        ax.set_xticks([])
        ax.set_xlabel(r"$x$")
        ax.set_yticks([])
        ax.set_ylabel(r"$y$")
        if ax_xz is not None:
            ax_xz.set_yticks([])
            ax_xz.set_ylabel(r"$z$")
        if ax_zy is not None:
            ax_zy.set_xticks([])
            ax_zy.set_xlabel(r"$z$")
    else:
        ax.set_xlabel(r"$x$ (%s)" % (po.r_unit_label))
        ax.set_ylabel(r"$y$ (%s)" % (po.r_unit_label))
        if ax_xz is not None:
            ax_xz.set_ylabel(r"$z$ (%s)" % (po.r_unit_label))
        if ax_zy is not None:
            ax_zy.set_xlabel(r"$z$ (%s)" % (po.r_unit_label))

    auto_nice_plot(ax, do_far_ticks=scale is None)
    if ax_xz is not None:
        auto_nice_plot(ax_xz, do_far_ticks=scale is None)
    if ax_zy is not None:
        auto_nice_plot(ax_zy, do_far_ticks=scale is None)

    # Save
    Fp_save = "%s/%s/%s_pos_trace_%d" % (po.dir_save, rsim.name, rsim.name, idx_orb)
    if po.Fp_suffix is not None:
        Fp_save += "%s" % po.Fp_suffix
    if ut.do_paper:
        savefig(Fp_save + ".pdf", dpi=400)
    else:
        savefig(Fp_save + ".png", dpi=po.dpi)
    print('\nSaved "%s" ' % Fp_save[-64:])
    plt.close()


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR rebound_plots.py  ====\n")
