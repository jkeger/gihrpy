"""
Functions for GIHR plots of snapshot data, etc.

Notation etc
------------
+ See README.md.
+ param_* is usually a parameter name (str), while p_* is a value (e.g. float).
"""

from jkeger import *
import gihrpy.utilities as ut
from gihrpy.classes_sim import Simulation, SimSet
from gihrpy.objects import Di_sim, Di_sim_set


def _plot_snapshot(
    A2_pos, A1_rgba, size=2**2, proj="xy", thick=None, A1_h=None, ax=None
):
    """Plot the particles from a snapshot, on existing axes.

    Parameters
    ----------
    A2_pos : [[float]]
        The particle positions [[x0, y0, z0], [x1, y1, z1], ...] (r_unit).

    A1_rgba : [[float]]
        The RGBA colour for each particle.

    size : float
        The marker size.

    proj : str
        The plotting projection plane: xy, xz, zy, or yz.

    thick : float
        The thickness of a slice in z about z=0 for plotting particles (r_unit).
        None    Plot all particles.
        -ve     Plot all particles between z_min and z=0 in |thick| slices.
        0       Plot the particles within their smoothing lengths of z=0.

    A1_h : float
        The particle smoothing lengths (r_unit), required for e.g. thick=0.

    ax : Axes
        The axes on which to plot.
    """
    if ax is None:
        ax = plt.gca()
    plt.sca(ax)

    if len(A2_pos) == 0:
        print("\n # Warning: no particles provided!")

    # Downsample to a subset of particles
    f_downsample = ut.extract_bonus("downsample", type=float)
    if f_downsample is not None:
        # Downsample to this fraction
        n_downsample = int(np.round(1 / f_downsample))
        # Keep every n^th particle
        A2_pos = A2_pos[::n_downsample]
        A1_rgba = A1_rgba[::n_downsample]
        if A1_h is not None:
            A1_h = A1_h[::n_downsample]

    # Projection axes
    if proj == "xy":
        x_co = 0
        y_co = 1
        z_co = 2
    elif proj == "xz":
        x_co = 0
        y_co = 2
        z_co = 1
    elif proj == "zy":
        x_co = 2
        y_co = 1
        z_co = 0
    elif proj == "yz":
        x_co = 1
        y_co = 2
        z_co = 0

    # Project xz view as if from -y not +y, for more intuitive zorder
    if proj == "xz":
        A2_pos[:, z_co] = -A2_pos[:, z_co]

    # Rasterise anything with zorder < 1
    ax.set_rasterization_zorder(1)

    # Always want to plot higher-z particles on top of lower ones
    A1_sort = A2_pos[:, z_co].argsort()
    A2_pos = A2_pos[A1_sort]
    A2_pos = A2_pos[::-1]
    A1_rgba = A1_rgba[A1_sort]
    A1_rgba = A1_rgba[::-1]

    # z slice(s)
    if thick is None:
        # Plot all particles, but still in slices in order of increasing z
        A1_z = np.linspace(np.amin(A2_pos[:, z_co]), np.amax(A2_pos[:, z_co]), 500)
    elif thick > 0:
        # Plot one slice below z=0
        A1_z = np.array([-thick, 0])
    elif thick < 0:
        # Plot slices of increasing z up to z=0
        z_min = np.amin(A2_pos[:, z_co])
        if z_min > 0:
            print(
                "\n # Warning: z_min (%.4g) > 0 for thick (%.4g) < 0" % (z_min, thick)
            )
            return

        A1_z = np.linspace(z_min, 0, int(np.floor(z_min / thick)))

        # Prevent a huge list of slices
        if 1000 < len(A1_z):
            A1_z = np.append(np.linspace(z_min, A1_z[-100], 900), A1_z[-99:])
    elif thick == 0:
        # Plot particles with |z| < h
        # Placeholder for a single slice
        A1_z = np.array([0, 0])

    # ========
    # Loop over each z slice
    # ========
    for i_z in range(len(A1_z) - 1):
        # Select particles by z
        if thick == 0:
            A1_sel_z = np.where(abs(A2_pos[:, z_co]) < A1_h)[0]
        else:
            A1_sel_z = np.where(
                (A1_z[i_z] < A2_pos[:, z_co]) & (A2_pos[:, z_co] < A1_z[i_z + 1])
            )[0]
        if len(A1_sel_z) == 0:
            continue

        A2_pos_z = A2_pos[A1_sel_z]
        A1_rgba_z = A1_rgba[A1_sel_z]

        # Convert colours to an array of tuples
        A1_rgba_z = A1_rgba_z.T
        A1_rgba_z = list(zip(A1_rgba_z[0], A1_rgba_z[1], A1_rgba_z[2], A1_rgba_z[3]))

        # Plot!
        ax.scatter(
            A2_pos_z[:, x_co],
            A2_pos_z[:, y_co],
            c=A1_rgba_z,
            edgecolors="none",
            marker=".",
            s=size,
            zorder=-2 - 2 * len(A1_z) + i_z,
        )


def plot_snapshots(
    sim, param_c, time=None, time_end=None, time_step=None, A1_fig_ax=None
):
    """Call _plot_snapshot() for one or more snapshots from a simulation.

    Parameters
    ----------
    sim : Simulation
        An object with information about a simulation.

    param_c : str
        How the particles are coloured. Most particle properties are valid
        (see Simulation.load_snapshot_data), plotted with colour maps or e.g. by
        material colour.

    time, time_end, time_step : int (opt.)
        The snapshot times to plot. See Simulation.select_snapshot_times.

    A1_fig_ax : [fig, ax] (opt.)
        If provided, then plot on this existing figure and axes instead of
        making new ones. Also assumes that only one snapshot is selected.

    Bonus options (set in ut.A1_bonus)
    -------------
    xz_panel        Add an x--z projection panel above the main x--y one.
    axes_off        Hide axes, e.g. to stop tick labels changing the size.
    fof_orbits      Plot the estimated orbits of fof groups.
    colour_via_fof  Colour particles by their fof group's mean/com param_c.
    init_orbits     Plot the initial orbits e.g. for a rebound collision.
    downsample=<f>  Downsample the particles, to a fraction f of the total.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim, Simulation):
        sim = Di_sim[sim]
    sim.init()
    po = sim.po

    time = check_none_or_def(time, None)
    time_end = check_none_or_def(time_end, None)
    time_step = check_none_or_def(time_step, None)

    # Track the runtime elapsed
    clock_start = datetime.datetime.now()

    # Select snapshots
    A1_time = sim.select_snapshot_times(time, time_end, time_step)

    # Colourbar or not
    if param_c in ["mat_id", "fof_id"]:
        do_cbar = False
    else:
        do_cbar = True

    # Misc prep
    if po.Di_init_orbit is not None:
        o_0 = sim.init_orbit()
    if "init_orbits" in ut.A1_bonus:
        o_0_1, o_0_2, o_0_c = sim.init_orbits()

    # ========
    # Load and plot each snapshot
    # ========
    for i_time, time in enumerate(A1_time):
        # Save filename
        Fp_save = "%s/%s/%s_snap_%s_%06d" % (
            po.dir_save,
            sim.name,
            sim.name,
            param_c,
            time,
        )
        if po.proj != "xy":
            Fp_save += "_%s" % po.proj
        if param_c == "fof_id":
            Fp_save += "_%s" % sim.link_len
        if "fof_orbits" in ut.A1_bonus:
            Fp_save += "_fof_orb"
        if po.Fp_suffix is not None:
            Fp_save += "%s" % po.Fp_suffix
        if ut.do_paper:
            Fp_save += ".pdf"
        else:
            Fp_save += ".png"
        # Skip if already plotted
        if ut.do_skip_exist:
            if os.path.isfile(Fp_save):
                continue

        # Print progress
        clock_cur = datetime.datetime.now()
        print(
            '%s Plotting snapshot "%d" (%d of %d) '
            % (clock_cur.strftime("%H:%M:%S"), time, i_time + 1, len(A1_time)),
        )

        # Figure and axes
        if A1_fig_ax is None:
            fig = plt.figure()
            # Custom or standard axes
            if param_c == "Q_v_inf":
                # Overcomplicated grid axes layout for two colour bars
                n_gs = 100
                gs = mpl.gridspec.GridSpec(nrows=n_gs, ncols=n_gs)
                x_div = n_gs - 12
                y_div = int(n_gs * 0.36)

                ax = plt.subplot(gs[:, :x_div], aspect="equal")
                cax = plt.subplot(gs[y_div : n_gs - 2, x_div : x_div + 3])
                cax_2 = plt.subplot(gs[2:y_div, x_div : x_div + 3])

                # Spacing
                gs.update(wspace=1.5, hspace=5)
                plt.sca(ax)
            else:
                ax = fig.add_subplot(111, aspect="equal")
                cax = None

            # Add z projection axes
            if "xz_panel" in ut.A1_bonus:
                divider = make_axes_locatable(ax)
                divider.set_aspect(True)
                if po.xz_size is not None:
                    size = po.xz_size
                else:
                    size = "50%"
                ax_xz = divider.append_axes("top", size=size, sharex=ax, pad=0)
                ax_xz.set_aspect("equal", adjustable="datalim")

                # Force main axes in x--y plane
                po.proj = "xy"
            else:
                ax_xz = None

            # Add colour bar axes
            if do_cbar and cax is None:
                # Separate axis for the colour bar to not distort the aspect ratio
                if "xz_panel" in ut.A1_bonus:
                    # Make room for new full-height axes on the right
                    fig.set_figwidth(11)
                    fig.subplots_adjust(right=0.8)
                    cax = fig.add_axes([0.813, 0.118, 0.025, 0.8487])  # [L, B, W, H]
                else:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
        else:
            if "xz_panel" in ut.A1_bonus:
                fig, ax, ax_xz = A1_fig_ax
            else:
                fig, ax = A1_fig_ax
            plt.figure(fig.number)

        # Background colour
        ax.set_facecolor(po.rgba_bg)

        # Load the required particle data
        A1_id, A1_mat_id, A1_m, A2_pos, A1_p_c = sim.load_snapshot_data(
            time, ["id", "mat_id", "m", "pos", param_c]
        )

        # Override the particle colour parameter by their fof group's value
        if "colour_via_fof" in ut.A1_bonus:
            A1_fof_id = sim.load_snapshot_data(time, "fof_id")
            A1_p_c_fof = sim.load_fof_data(time, param_c)
            for fof_id in range(po.fof_id_max + 1):
                # Set all particles in this group to the group's value
                A1_sel_fof = np.where(A1_fof_id == fof_id)[0]
                A1_p_c[A1_sel_fof] = A1_p_c_fof[fof_id]

        # Convert to plotting units
        A2_pos /= po.r_unit
        if param_c == "Q_v_inf":
            # Q (po.r_unit) for bound particles, v_inf (km/s) for unbound
            oset = sim.load_snapshot_data(time, "o")
            A1_sel_bnd = np.array(oset.A1_e < 1)
            A1_sel_unb = np.array(oset.A1_e >= 1)
            A1_p_c[A1_sel_bnd] /= po.r_unit
            A1_p_c[A1_sel_unb] /= 1e3
        else:
            A1_p_c = A1_p_c / po.Di_param_unit_label[param_c].unit

        # Axis limits and ticks
        if po.func_ax_lim is not None:
            A1_ax_lim, tick_base, tick_base_minor = po.func_ax_lim(po, time, A2_pos)
        else:
            A1_ax_lim = po.A1_ax_lim
            tick_base = po.tick_base
            tick_base_minor = po.tick_base_minor

        # Ignore particles outside the axis limits, with a little padding
        dx_lim = A1_ax_lim[1] - A1_ax_lim[0]
        dy_lim = A1_ax_lim[3] - A1_ax_lim[2]
        pad = 0.02
        A1_sel_lim = np.array(
            (A2_pos[:, 0] > A1_ax_lim[0] - pad * dx_lim)
            & (A2_pos[:, 0] < A1_ax_lim[1] + pad * dx_lim)
            & (A2_pos[:, 1] > A1_ax_lim[2] - pad * dy_lim)
            & (A2_pos[:, 1] < A1_ax_lim[3] + pad * dy_lim)
        )

        # Scale marker size with axis size
        if po.func_marker_size is not None:
            marker_size = po.func_marker_size(po, A1_ax_lim)
        else:
            marker_size = po.marker_size

        # Colour map and bar
        if param_c == "Q_v_inf":
            po.Di_cmap_param["A1_sel_unb"] = A1_sel_unb[A1_sel_lim]
            (
                cmap,
                vmin,
                vmax,
                norm,
                extend,
                cmap_2,
                vmin_2,
                vmax_2,
                norm_2,
                extend_2,
            ) = po.func_cmap_param(po, param_c, A1_p_c[A1_sel_lim])

            # Plot the colour bar, unless plotting on existing axes
            if A1_fig_ax is None:
                scat = ax.scatter(
                    [],
                    [],
                    c=[],
                    cmap=cmap,
                    edgecolors="none",
                    marker=".",
                    vmin=vmin,
                    vmax=vmax,
                    norm=norm,
                    zorder=0,
                )
                scat_2 = ax.scatter(
                    [],
                    [],
                    c=[],
                    cmap=cmap_2,
                    edgecolors="none",
                    marker=".",
                    vmin=vmin_2,
                    vmax=vmax_2,
                    norm=norm_2,
                    zorder=0,
                )

                cbar = plt.colorbar(scat, cax=cax, extend=extend)
                cbar_2 = plt.colorbar(scat_2, cax=cax_2, extend=extend_2)

                cbar.set_label(sim.param_label_unit("Q"), size=20)
                cbar_2.set_label(r"Unbound $v_{\infty}$ (km~s$^{-1}$)", size=20)
                set_large_ticks(cax)
                set_large_ticks(cax_2)
                if cax_2.get_ylim()[1] > 0.4 and cax_2.get_ylim()[1] < 1:
                    cax_2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.2))
            else:
                cbar = None

            # Store
            Di_cmap = {
                "cmap": cmap,
                "norm": norm,
                "A1_sel": A1_sel_bnd,
                "cmap_2": cmap_2,
                "norm_2": norm_2,
                "A1_sel_2": A1_sel_unb,
            }
        elif do_cbar:
            cmap, vmin, vmax, norm, extend = po.func_cmap_param(
                po, param_c, A1_p_c[A1_sel_lim]
            )

            # Plot the colour bar, unless plotting on existing axes
            if A1_fig_ax is None:
                scat = ax.scatter(
                    [],
                    [],
                    c=[],
                    cmap=cmap,
                    edgecolors="none",
                    marker=".",
                    vmin=vmin,
                    vmax=vmax,
                    norm=norm,
                    zorder=0,
                )

                cbar = plt.colorbar(scat, cax=cax, extend=extend)

                cbar.set_label(sim.param_label_unit(param_c), size=po.fontsize)
                set_large_ticks(cax)
                if norm is None:
                    set_std_form_cbar(cbar, A1_p_c[A1_sel_lim])
            else:
                cbar = None

            # Store
            Di_cmap = {"cmap": cmap, "norm": norm}
        else:
            # Placeholder
            cbar = None
            Di_cmap = None

        # Set RGBA colours
        A1_rgba = po.colours_from_data(param_c, A1_p_c, Di_cmap)

        # Also load smoothing lengths if needed
        if po.thick == 0 or po.thick_xz == 0:
            A1_h = sim.load_snapshot_data(time, "h")

        # Plot!
        _plot_snapshot(
            A2_pos[A1_sel_lim],
            A1_rgba[A1_sel_lim],
            size=marker_size,
            proj=po.proj,
            thick=po.thick,
            A1_h=A1_h[A1_sel_lim] if po.thick == 0 else None,
            ax=ax,
        )

        # Axes etc
        ax.set_xlabel(r"$%s$ (%s)" % (po.proj[0], po.r_unit_label))
        ax.set_ylabel(r"$%s$ (%s)" % (po.proj[1], po.r_unit_label))
        ax.set_xlim(A1_ax_lim[0], A1_ax_lim[1])
        ax.set_ylim(A1_ax_lim[2], A1_ax_lim[3])
        if tick_base is not None:
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=tick_base))
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=tick_base))
        if tick_base_minor is not None:
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=tick_base_minor))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=tick_base_minor))

        # Turn off axes
        if "axes_off" in ut.A1_bonus:
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)

        # Plot x--z projection panel
        if "xz_panel" in ut.A1_bonus:
            # Axis limits
            if po.xz_size is not None:
                size_rel = float(po.xz_size[:-1]) / 100
            else:
                size_rel = 0.5
            z_lim = (A1_ax_lim[1] - A1_ax_lim[0]) / 2 * size_rel
            A1_ax_lim_xz = np.array([A1_ax_lim[0], A1_ax_lim[1], -z_lim, z_lim])

            # Ignore particles outside the axis limits, with a little padding
            A1_sel_lim_xz = np.array(
                (A2_pos[:, 0] > A1_ax_lim_xz[0] - pad * dx_lim)
                & (A2_pos[:, 0] < A1_ax_lim_xz[1] + pad * dx_lim)
                & (A2_pos[:, 1] > A1_ax_lim_xz[2] - pad * dy_lim)
                & (A2_pos[:, 1] < A1_ax_lim_xz[3] + pad * dy_lim)
            )

            # Plot!
            _plot_snapshot(
                A2_pos[A1_sel_lim_xz],
                A1_rgba[A1_sel_lim_xz],
                size=marker_size,
                proj="xz",
                thick=po.thick_xz,
                A1_h=A1_h[A1_sel_lim_xz] if po.thick_xz == 0 else None,
                ax=ax_xz,
            )

            # Axes etc
            ax_xz.set_ylabel(sim.param_label_unit("z"))
            ax_xz.set_ylim(A1_ax_lim_xz[2], A1_ax_lim_xz[3])
            if tick_base is not None:
                ax_xz.yaxis.set_major_locator(
                    mpl.ticker.MultipleLocator(base=tick_base)
                )
            if tick_base_minor is not None:
                ax_xz.yaxis.set_minor_locator(
                    mpl.ticker.MultipleLocator(base=tick_base_minor)
                )
            plt.setp(ax_xz.get_xticklabels(), visible=False)
            plt.setp(ax_xz.get_xticklines(), visible=False)
            ax_xz.axhline(0, c="k", lw=0.3, zorder=-999999)

            # Turn off axes
            if "axes_off" in ut.A1_bonus:
                ax_xz.axes.yaxis.set_visible(False)

            # Misc extras
            if po.Di_init_orbit is not None:
                po.plot_init_orbit(o_0, ax=ax_xz)
            if "xp" in sim.impact.A1_preset_1 and "planet" in po.Di_misc_colour.keys():
                sim.plot_planet(A1_ax_lim=A1_sel_lim_xz, po=po, ax=ax_xz)

            auto_nice_plot(ax, cbar=None, do_far_ticks="axes_off" not in ut.A1_bonus)

        # Add text annotations etc
        sim.add_text(time, ax=ax)

        # Misc extras
        if po.Di_init_orbit is not None:
            po.plot_init_orbit(o_0, ax=ax)
        if "xp" in sim.impact.A1_preset_1 and "planet" in po.Di_misc_colour.keys():
            sim.plot_planet(A1_ax_lim=A1_ax_lim, po=po, ax=ax)
        if po.func_snapshot_add_misc is not None:
            po.func_snapshot_add_misc(po=po, time=time, param_c=param_c, ax=ax)

        # Inset panels
        A1_ax_in = []
        for po_in in po.A1_po_inset:
            # Inset axes
            ax_in = inset_axes(ax, width=0, height=0)
            ip = InsetPosition(ax, po_in.A1_inset_loc)
            ax_in.set_axes_locator(ip)
            A1_ax_in.append(ax_in)

            # Background colour, set slightly transparent
            rgba_bg = list(po_in.rgba_bg)
            rgba_bg[3] *= 0.85
            ax_in.set_facecolor(tuple(rgba_bg))

            # Axis limits and ticks
            if po_in.func_ax_lim is not None:
                A1_ax_lim_in, tick_base_in, tick_base_minor_in = po_in.func_ax_lim(
                    po_in, time, A2_pos
                )
            else:
                A1_ax_lim_in = po_in.A1_ax_lim
                tick_base_in = po_in.tick_base
                tick_base_minor_in = po_in.tick_base_minor
            ax_in.set_xlim(A1_ax_lim_in[0], A1_ax_lim_in[1])
            ax_in.set_ylim(A1_ax_lim_in[2], A1_ax_lim_in[3])

            # Ignore particles outside the axis limits, with a little padding
            dx_lim = A1_ax_lim_in[1] - A1_ax_lim_in[0]
            dy_lim = A1_ax_lim_in[3] - A1_ax_lim_in[2]
            pad = 0.02
            A1_sel_lim_in = np.array(
                (A2_pos[:, 0] > A1_ax_lim_in[0] - pad * dx_lim)
                & (A2_pos[:, 0] < A1_ax_lim_in[1] + pad * dx_lim)
                & (A2_pos[:, 1] > A1_ax_lim_in[2] - pad * dy_lim)
                & (A2_pos[:, 1] < A1_ax_lim_in[3] + pad * dy_lim)
            )

            # Scale marker size with axis size
            if po_in.func_marker_size is not None:
                marker_size_in = po_in.func_marker_size(po_in, A1_ax_lim_in)
            else:
                marker_size_in = po_in.marker_size

            # Also load smoothing lengths if needed
            if po_in.thick == 0:
                A1_h = sim.load_snapshot_data(time, "h")

            # Plot!
            _plot_snapshot(
                A2_pos[A1_sel_lim_in],
                A1_rgba[A1_sel_lim_in],
                size=marker_size_in,
                proj=po_in.proj,
                thick=po_in.thick,
                A1_h=A1_h[A1_sel_lim_in] if po_in.thick == 0 else None,
                ax=ax_in,
            )

            # Axes, no ticks or labels
            ax_in.axes.xaxis.set_visible(False)
            ax_in.axes.yaxis.set_visible(False)
            ax_in.set_xlim(A1_ax_lim_in[0], A1_ax_lim_in[1])
            ax_in.set_ylim(A1_ax_lim_in[2], A1_ax_lim_in[3])
            for side in ["bottom", "top", "left", "right"]:
                ax_in.spines[side].set_color("0.4")
                ax_in.spines[side].set_linewidth(1)

            # Mark inset
            po_in.mark_inset_panel(ax, ax_in)

            # Misc extras
            if po_in.Di_init_orbit is not None:
                po_in.plot_init_orbit(o_0, ax=ax_in)
            if (
                "xp" in sim.impact.A1_preset_1
                and "planet" in po_in.Di_misc_colour.keys()
            ):
                sim.plot_planet(A1_ax_lim=A1_ax_lim_in, po=po_in, ax=ax_in)
            if po_in.func_snapshot_add_misc is not None:
                po_in.func_snapshot_add_misc(
                    po=po_in, time=time, param_c=param_c, ax=ax
                )

        # Bonus options
        if "fof_orbits" in ut.A1_bonus:
            oset_fof = sim.load_or_compute_oset_fof(time, sim.link_len, po.fof_id_max)

            # Colour by ID or parameter
            if "colour_via_fof" in ut.A1_bonus:
                A1_c_fof = np.empty((po.fof_id_max + 1, 4))
                A1_fof_id = sim.load_snapshot_data(time, "fof_id")
                for fof_id in range(po.fof_id_max + 1):
                    # Extract the already-set colour from any particle in this group
                    A1_sel_fof = np.where(A1_fof_id == fof_id)[0]
                    A1_c_fof[fof_id] = A1_rgba[A1_sel_fof[0]]
            else:
                A1_c_fof = po.cmap(np.linspace(0, 1, po.fof_id_max + 1))

            # Plot each fof group
            for fof_id, (o, colour) in enumerate(zip(oset_fof.A1_o, A1_c_fof)):
                # Set zorder, rasterise if many orbits
                if po.fof_id_max > 16:
                    zorder = -0.1 + (po.fof_id_max - fof_id) / po.fof_id_max
                else:
                    zorder = 2 * po.fof_id_max - fof_id

                # Bound or unbound linestyle
                if o.e < 1:
                    ls = "-"
                else:
                    ls = "--"

                # Plot the orbit
                for ax_i, proj in [[ax, po.proj], [ax_xz, "xz"]]:
                    if ax_i is None:
                        continue

                    plot_orbit(
                        o,
                        A1_pos_0=sim.A1_pos_p,
                        unit=po.r_unit,
                        proj=proj,
                        gradient="nu",
                        c=colour,
                        lw=0.8,
                        alpha=0.8,
                        ls=ls,
                        A1_ax_lim="auto",
                        mark_periapsis=False,
                        zorder=zorder,
                        ax=ax_i,
                    )
        if "init_orbits" in ut.A1_bonus:
            # Plot initial and centre-of-mass orbits
            for o_i, c_i in [
                [o_0_1, po.Di_misc_colour["o_1"]],
                [o_0_2, po.Di_misc_colour["o_2"]],
                [o_0_c, po.Di_misc_colour["com"]],
            ]:
                for ax_i, proj in [[ax, po.proj], [ax_xz, "xz"]]:
                    if ax_i is None:
                        continue

                    plot_orbit(
                        o_i,
                        A1_pos_0=[0, 0, 0],
                        unit=po.r_unit,
                        proj=proj,
                        c=c_i,
                        lw=0.5,
                        alpha=0.6,
                        ls="-",
                        num_step=int(1e6),
                        A1_ax_lim="auto",
                        zorder=-999,
                        ax=ax_i,
                    )

        # Ensure main axis limits haven't changed
        ax.set_xlim(A1_ax_lim[0], A1_ax_lim[1])
        ax.set_ylim(A1_ax_lim[2], A1_ax_lim[3])

        auto_nice_plot(
            ax,
            cbar=cbar,
            fontsize=po.fontsize,
            do_far_ticks="axes_off" not in ut.A1_bonus,
        )

        # Save and close, unless adding to an existing figure
        if A1_fig_ax is None:
            savefig(Fp_save, dpi=po.dpi)
            plt.close()
        else:
            print("")
            return

    print('\nSaved "%s" ' % Fp_save[-64:])
    plt.close()

    # Print the runtime elapsed
    clock_end = datetime.datetime.now()
    print(
        "Finished at %s, elapsed %s"
        % (clock_end.strftime("%H:%M:%S"), str(clock_end - clock_start).split(".")[0])
    )


def plot_particle_params(
    sim_set, param_y, param_x, param_c=None, time_all=None, A1_fig_ax=None
):
    """Plot one particle parameter against another.

    From one or more simulations and snapshots.

    e.g. radial profile, phase diagram, orbital elements.

    Parameters
    ----------
    sim_set : Simulation or SimSet
        An object with information about one or a set of simulations.

    param_y, param_x : str
        The particle properties to plot on the y and x axes.

    param_c : str (opt.)
        How the particles are coloured. Most particle properties are valid
        (see Simulation.load_snapshot_data), plotted with colour maps or e.g. by
        material colour. Or None for a single colour (per simulation).

    time_all : int
        The time of the snapshot to plot (s), if not already set by sim_set.
        0       The initial snapshot.
        -1      The final snapshot.

    A1_fig_ax : [fig, ax]
        If provided, then plot on this existing figure and axes instead of
        making new ones.

    Bonus options (set in ut.A1_bonus)
    -------------
    init_prof       Plot the initial profile(s) as well as the particles.
    init_orbits     Plot the initial orbits e.g. for a rebound collision.
    downsample=<f>  Downsample the particles, to a fraction f of the total.
    fof_only        Instead of individual particles, plot fof group parameters.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim_set, Simulation) and not isinstance(sim_set, SimSet):
        try:
            sim_set = Di_sim[sim_set]
        except KeyError:
            sim_set = Di_sim_set[sim_set]
    # Wrap a single Simulation as a SimSet object if needed
    if isinstance(sim_set, Simulation):
        is_single_sim = True
        sim_set = SimSet(
            "%s" % sim_set.name,
            A1_sim=[sim_set],
            A1_colour=[A1_c[0]],
            A1_linestyle=["-"],
            A1_label=[None],
        )
    else:
        is_single_sim = False
    sim_set.init()

    param_c = check_none_or_def(param_c, None)
    time_all = check_none_or_def(time_all, None)
    if time_all is not None:
        sim_set.A1_time = [time_all] * sim_set.num_sim

    # Colourbar or not
    if param_c in [None, "mat_id"]:
        do_cbar = False
    else:
        do_cbar = True

    # Figure
    if A1_fig_ax is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        fig, ax = A1_fig_ax
        plt.figure(fig.number)

    # ========
    # Load and plot each simulation
    # ========
    for i_sim, sim in enumerate(sim_set.A1_sim):
        sim.init()
        po = sim.po

        # Select snapshot
        time = sim.select_snapshot_times(sim_set.A1_time[i_sim])[0]

        # Load the required particle (or fof) data
        if "fof_only" in ut.A1_bonus:
            A1_p_x, A1_p_y, A1_p_c = sim.load_fof_data(
                time, [param_x, param_y, param_c]
            )
        else:
            A1_p_x, A1_p_y, A1_p_c = sim.load_snapshot_data(
                time, [param_x, param_y, param_c]
            )

        # Plotting units
        A1_p_x = A1_p_x / po.Di_param_unit_label[param_x].unit
        A1_p_y = A1_p_y / po.Di_param_unit_label[param_y].unit
        A1_p_c = A1_p_c / po.Di_param_unit_label[param_c].unit

        # Restrict and sort by x
        if param_x in po.Di_param_lim.keys():
            x_lim_0, x_lim_1 = po.Di_param_lim[param_x]
            x_pad = 0.02 * (x_lim_1 - x_lim_0)
            A1_sel_x = np.where(
                (x_lim_0 - x_pad < A1_p_x) & (A1_p_x < x_lim_1 + x_pad)
            )[0]
            A1_p_x = A1_p_x[A1_sel_x]
            A1_p_y = A1_p_y[A1_sel_x]
            A1_p_c = A1_p_c[A1_sel_x]
            if len(A1_p_x) == 0:
                print("\n # Warning: no particles selected")
        else:
            x_lim_0 = np.nanmin(A1_p_x)
            x_lim_1 = np.nanmax(A1_p_x)
        A1_sort_x = np.argsort(A1_p_x)
        A1_p_x = A1_p_x[A1_sort_x]
        A1_p_y = A1_p_y[A1_sort_x]
        A1_p_c = A1_p_c[A1_sort_x]

        # Restrict to a subset of particles
        f_downsample = ut.extract_bonus("downsample", type=float)
        if f_downsample is not None:
            # Downsample to this fraction
            n_downsample = int(np.round(1 / f_downsample))
            # Keep every n^th particle
            A1_p_x = A1_p_x[::n_downsample]
            A1_p_y = A1_p_y[::n_downsample]
            A1_p_c = A1_p_c[::n_downsample]
        if param_x == "r_xy" and param_y != "period":
            # Restrict to thin equatorial slice
            A1_sel_z = np.where(abs(A2_pos[:, 2]) < A1_h)[0]
            A1_p_x = A1_p_x[A1_sel_z]
            A1_p_y = A1_p_y[A1_sel_z]
            A1_p_c = A1_p_c[A1_sel_z]

        # Colour bar
        cmap, vmin, vmax, norm, extend = po.func_cmap_param(po, param_c, A1_p_c)
        Di_cmap = {"cmap": cmap, "norm": norm}
        if cmap is not None:
            scat = ax.scatter(
                [],
                [],
                c=[],
                cmap=cmap,
                edgecolors="none",
                marker=".",
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                zorder=0,
            )

            cbar = plt.colorbar(scat, extend=extend)
            cbar.set_label(sim.param_label_unit(param_c), size=po.fontsize)
            set_large_ticks(cbar.ax)
            if norm is None:
                set_std_form_cbar(cbar, A1_p_c)
        else:
            # Placeholder
            cbar = None

        # Set colours
        if param_c is not None:
            A1_colour = po.colours_from_data(param_c, A1_p_c, Di_cmap)
        else:
            # Single colour (per simulation)
            A1_colour = sim_set.A1_colour[i_sim]

        # Rasterise anything with zorder < 1
        ax.set_rasterization_zorder(1)

        # Scatter properties
        alpha = 0.2
        size = 2**2
        if sim_set.A1_zorder is not None:
            zorder = -2 + sim_set.A1_zorder[i_sim] / np.amax(sim_set.A1_zorder)
        else:
            zorder = -2

        # Plot!
        plt.scatter(
            A1_p_x,
            A1_p_y,
            c=A1_colour,
            edgecolor="none",
            marker=".",
            s=size,
            alpha=alpha,
            zorder=zorder,
        )

        # Misc extras
        if "init_prof" in ut.A1_bonus:
            ut.plot_particle_params_init_prof(sim, param_x, param_y, ax)
        if "init_orbits" in ut.A1_bonus:
            ut.plot_particle_params_init_orbits(sim, param_x, param_y, ax)
        if po.func_particle_params_extra is not None:
            po.func_particle_params_extra(sim, time, param_x, param_y, ax)

        # Legend entries
        if not is_single_sim and sim_set.A1_label is not None:
            plt.scatter(
                [],
                [],
                c=sim_set.A1_colour[i_sim],
                marker=".",
                s=13**2,
                label=sim_set.A1_label[i_sim],
            )

    # Legend
    if not is_single_sim and sim_set.A1_label is not None:
        leg = plt.legend(prop={"size": po.fontsize_text}, title=sim_set.legend_title)
        plt.setp(leg.get_title(), fontsize=po.fontsize)

    # Reference simulation for label utilities etc
    sim = sim_set.A1_sim[0]
    po = sim_set.po

    # ========
    # Axes etc
    # ========
    ax.set_xlim(x_lim_0, x_lim_1)
    if param_y in po.Di_param_lim.keys():
        y_lim_0, y_lim_1 = po.Di_param_lim[param_y]
        ax.set_ylim(y_lim_0, y_lim_1)
    else:
        set_auto_limits(ax, A1_p_x, A1_p_y, set_x=False, set_y=True)

    ax.set_xlabel(sim.param_label_unit(param_x))
    ax.set_ylabel(sim.param_label_unit(param_y))

    if param_x in po.Di_param_log.keys() and po.Di_param_log[param_x]:
        ax.set_xscale("log")
    if param_y in po.Di_param_log.keys() and po.Di_param_log[param_y]:
        ax.set_yscale("log")

    auto_nice_plot(ax, cbar=cbar, fontsize=po.fontsize)

    # Save and close, unless adding to an existing figure
    if A1_fig_ax is None:
        if is_single_sim:
            Fp_save = "%s/%s/%s" % (po.dir_save, sim_set.name, sim_set.name)
        else:
            Fp_save = "%s/%s" % (po.dir_save, sim_set.name)
        Fp_save += "_pp_%s_%s_%s" % (
            param_y,
            param_x,
            param_c,
        )
        if time_all is not None:
            Fp_save += "_%06d" % sim.select_snapshot_times(time_all)[0]
        if param_c == "fof_id":
            Fp_save += "_%s" % sim.link_len
        if po.Fp_suffix is not None:
            Fp_save += "%s" % po.Fp_suffix

        # Save
        if ut.do_paper:
            savefig(Fp_save + ".pdf", dpi=400)
        else:
            savefig(Fp_save + ".png", dpi=po.dpi)
        print('\nSaved "%s" ' % Fp_save[-64:])

        plt.close()
    else:
        print("")
        return


def plot_param_hist_cum(
    sim_set, param_x, hist_type=None, time_all=None, do_cum_gt_x=False, A1_fig_ax=None
):
    """Plot a cumulative histogram of a parameter (with no binning).

    From one or more simulations and snapshots.

    Parameters
    ----------
    sim_set : Simulation or SimSet
        An object with information about one or a set of simulations.

    param_x : str
        The particle property for the cumulative histogram along the x axis.

    hist_type : str
        The type of histogram to set the y axis, choose from:
        frac        The cumulative fraction of the total number (default).
        count       The cumulative number.
        mass        The cumulative mass.
        mass_frac   The cumulative fraction of the total mass.

    do_cum_gt_x : bool
        If True, then plot the cumulative data above x rather than below.

    time_all : int
        The time of the snapshot to plot (s), if not already set by sim_set.
        0       The initial snapshot.
        -1      The final snapshot.

    A1_fig_ax : [fig, ax]
        If provided, then plot on this existing figure and axes instead of
        making new ones.

    Bonus options (set in ut.A1_bonus)
    -------------
    fof_only        Instead of individual particles, plot fof groups.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim_set, Simulation) and not isinstance(sim_set, SimSet):
        try:
            sim_set = Di_sim[sim_set]
        except KeyError:
            sim_set = Di_sim_set[sim_set]
    # Wrap a single Simulation as a SimSet object if needed
    if isinstance(sim_set, Simulation):
        is_single_sim = True
        sim_set = SimSet(
            "%s" % sim_set.name,
            A1_sim=[sim_set],
            A1_colour=[A1_c[0]],
            A1_linestyle=["-"],
            A1_label=[None],
        )
    else:
        is_single_sim = False
    sim_set.init()

    hist_type = check_option_or_def(
        hist_type, ["frac", "count", "mass", "mass_frac"], "frac"
    )
    time_all = check_none_or_def(time_all, None)
    if time_all is not None:
        sim_set.A1_time = [time_all] * sim_set.num_sim
    do_cum_gt_x = check_bool_or_def(do_cum_gt_x, False)

    # Figure
    if A1_fig_ax is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        fig, ax = A1_fig_ax
        plt.figure(fig.number)

    # ========
    # Load and plot each simulation
    # ========
    for i_sim, sim in enumerate(sim_set.A1_sim):
        sim.init()
        po = sim.po

        # Select snapshot
        time = sim.select_snapshot_times(sim_set.A1_time[i_sim])[0]

        # Load the required particle (or fof) data
        if hist_type in ["mass", "mass_frac"]:
            if "fof_only" in ut.A1_bonus:
                A1_p_x, A1_m = sim.load_fof_data(time, [param_x, "m"])
            else:
                A1_p_x, A1_m = sim.load_snapshot_data(time, [param_x, "m"])
        else:
            if "fof_only" in ut.A1_bonus:
                A1_p_x = sim.load_fof_data(time, param_x)
            else:
                A1_p_x = sim.load_snapshot_data(time, param_x)

        # Plotting units
        A1_p_x = A1_p_x / po.Di_param_unit_label[param_x].unit

        # Sort
        A1_sort_x = np.argsort(A1_p_x)
        if do_cum_gt_x:
            # Reverse order
            A1_sort_x = A1_sort_x[::-1]

        A1_p_x = A1_p_x[A1_sort_x]
        if hist_type in ["mass", "mass_frac"]:
            A1_m = A1_m[A1_sort_x]

        # Cumulative data
        if hist_type in ["frac", "count"]:
            A1_p_cum = np.arange(len(A1_p_x)) + 1
        elif hist_type in ["mass", "mass_frac"]:
            A1_p_cum = np.cumsum(A1_m)

        # Normalise
        if hist_type == "frac":
            A1_p_cum /= len(A1_p_x)
        elif hist_type == "mass_frac":
            A1_p_cum /= sum(A1_m)

        # Limits
        x_lim_0 = 0.9 * np.amin(A1_p_x)
        x_lim_1 = 1.1 * np.amax(A1_p_x)

        # Extend the start and end lines for clarity
        A1_x_step = np.append(
            A1_p_x[0], np.append(A1_p_x, x_lim_0 if do_cum_gt_x else x_lim_1)
        )
        A1_y_step = np.append(0, np.append(A1_p_cum, A1_p_cum[-1]))

        # Plot!
        ax.step(
            A1_x_step,
            A1_y_step,
            where="pre" if do_cum_gt_x else "post",
            c=sim_set.A1_colour[i_sim],
            lw=1.7,
            ls=sim_set.A1_linestyle[i_sim],
            alpha=0.8,
        )

        # Legend entries
        if not is_single_sim and sim_set.A1_label is not None:
            plt.scatter(
                [],
                [],
                c=sim_set.A1_colour[i_sim],
                lw=1.7,
                ls=sim_set.A1_linestyle[i_sim],
                label=sim_set.A1_label[i_sim],
            )

    # Reference simulation for label utilities etc
    sim = sim_set.A1_sim[0]
    po = sim_set.po

    # Legend
    if not is_single_sim and sim_set.A1_label is not None:
        leg = plt.legend(prop={"size": po.fontsize_text}, title=sim_set.legend_title)
        plt.setp(leg.get_title(), fontsize=po.fontsize)

    # ========
    # Axes etc
    # ========
    ax.set_xlim(x_lim_0, x_lim_1)
    if do_cum_gt_x:
        x_unit_label = po.Di_param_unit_label[param_x]
        x_sym = x_unit_label.symbol
        ax.set_xlabel(
            r"%s, %s (%s)" % (x_unit_label.label, x_sym, x_unit_label.unit_label)
        )
    else:
        ax.set_xlabel(sim.param_label_unit(param_x))
    if param_x in po.Di_param_log.keys() and po.Di_param_log[param_x]:
        ax.set_xscale("log")

    if hist_type == "frac":
        ax.set_ylim(0, 1)
        if do_cum_gt_x:
            ax.set_ylabel(r"Cumulative count, $f(>\!%s)$" % x_sym)
        else:
            ax.set_ylabel(r"Cumulative fraction")
    elif hist_type == "count":
        ax.set_ylim(0, None)
        if do_cum_gt_x:
            ax.set_ylabel(r"Cumulative count, $N(>\!%s)$" % x_sym)
        else:
            ax.set_ylabel(r"Cumulative count")
        if A1_p_cum[-1] > 300:
            ax.set_yscale("log")
    elif hist_type == "mass":
        ax.set_ylim(0, None)
        if do_cum_gt_x:
            ax.set_ylabel(r"Cumulative count, $m(>\!%s)$" % x_sym)
        else:
            ax.set_ylabel(r"Cumulative mass (kg)")
        if "m" in po.Di_param_log.keys() and po.Di_param_log["m"]:
            ax.set_yscale("log")
    elif hist_type == "mass_frac":
        ax.set_ylim(0, 1)
        if do_cum_gt_x:
            ax.set_ylabel(r"Cumulative count, $f_m(>\!%s)$" % x_sym)
        else:
            ax.set_ylabel(r"Cumulative mass fraction")

    auto_nice_plot(ax, fontsize=po.fontsize)

    # Save and close, unless adding to an existing figure
    if A1_fig_ax is None:
        if is_single_sim:
            Fp_save = "%s/%s/%s" % (po.dir_save, sim_set.name, sim_set.name)
        else:
            Fp_save = "%s/%s" % (po.dir_save, sim_set.name)
        Fp_save += "_hist_cum_%s_%s" % (hist_type, param_x)
        if time_all is not None:
            Fp_save += "_%06d" % sim.select_snapshot_times(time_all)[0]
        if po.Fp_suffix is not None:
            Fp_save += "%s" % po.Fp_suffix

        # Save
        if ut.do_paper:
            savefig(Fp_save + ".pdf", dpi=400)
        else:
            savefig(Fp_save + ".png", dpi=po.dpi)
        print('\nSaved "%s" ' % Fp_save[-64:])

        plt.close()
    else:
        print("")
        return


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR snapshot_plots.py  ====\n")
