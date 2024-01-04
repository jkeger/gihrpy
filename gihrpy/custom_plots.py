"""
Functions for GIHR customisable plots, etc.
"""

from jkeger import *
import gihrpy.utilities as ut
from gihrpy.classes_sim import Simulation, SimSet
from gihrpy.classes_reb import RebSim
from gihrpy.objects import Di_sim_set
from gihrpy.init_cond import plot_init_prof
from gihrpy.snapshot_plots import (
    plot_snapshots,
    plot_particle_params,
    plot_param_hist_cum,
)
from gihrpy.scenario_plots import plot_sim_param_subsets
from gihrpy.rebound_plots import reb_plot_orbit_evol


def plot_custom_tile(
    num_x,
    num_y,
    Fp_save,
    *args,
    x_space=None,
    y_space=None,
    share_x_l=True,
    share_y_l=True,
    share_x_tl=True,
    share_y_tl=True,
):
    """Make a tiled figure of panels using other plotting functions.

    Awkward way to do the separate function arguments but it works...

    Parameters
    ----------
    num_x, num_y : int
        The numbers of horizontal and vertical panels.

    Fp_save : str
        The filename (full path) for saving the figure.

    *args : ! ... ! ... !
        The arguments for each panel, preceeded, separated, and followed by "!"
        symbols. i.e. the same text that would be used on the command line with
        gihr.py, starting with the function to call. Order from the top row,
        left to right, then the 2nd row, etc. The functions to be called must
        have A1_fig_ax as the final argument.

    x_space, y_space : float (opt.)
        The x and y spacing between the panels, can be set negative.

    share_x_l, share_y_l, share_x_tl, share_y_tl : bool (opt.)
        If True then remove the x/y labels and/or x/y tick labels, respectively,
        apart from left- and bottom-edge panels.

    Bonus options (set in ut.A1_bonus)
    -------------
    no_min_x_tl, no_min_y_tl
        Remove the lowest x or y tick label, e.g. to help fit with another row
        of panels underneath this one.

    hide_shared_x_t, hide_shared_y_t
        Remove ticks on shared axes as well as tick labels.

    figsize=[float,float]
        An override figsize.

    dpi=int
        The dots-per-inch resolution. May want a lower value with many panels.

    panel_labels=*
        Overlay text labels for each panel in the specified location:
        "<ul><lcr><aA1iI>" for the upper/lower left/centre/right location and
        a,b,c,/A,B,C,/1,2,3,/i,ii,iii,/I,II,III,... labels.

    ix_match_y_lim=*, iy_match_x_lim=*
        Copy the y or x axis limits from this index axis, if share_y,x_tl is
        True. E.g. ix_match_y_lim=0 to copy the leftmost y limits.
    """
    num_x = int(num_x)
    num_y = int(num_y)
    num_tile = num_x * num_y
    A1_arg = np.asarray(args)
    ut.do_paper = True

    # Separate the function arguments for each panel
    A2_arg = [[] for i in range(num_tile)]
    A1_idx = np.where(A1_arg == "!")[0]
    assert (
        len(A1_idx) == num_tile + 1
    ), 'Must give num_tile=%d sets of arguments between %d "!" symbols' % (
        num_tile,
        num_tile + 1,
    )
    for i in range(num_tile):
        A2_arg[i] = list(A1_arg[A1_idx[i] + 1 : A1_idx[i + 1]])

    # Parse remaining (optional) arguments
    A1_opt = A1_arg[A1_idx[-1] + 1 :]
    num_opt = len(A1_opt)
    if num_opt > 0:
        x_space = check_float_or_def(A1_opt[0], None)
    if num_opt > 1:
        y_space = check_float_or_def(A1_opt[1], None)
    if num_opt > 2:
        share_x_l = check_bool_or_def(A1_opt[2], True)
    if num_opt > 3:
        share_y_l = check_bool_or_def(A1_opt[3], True)
    if num_opt > 4:
        share_x_tl = check_bool_or_def(A1_opt[4], True)
    if num_opt > 5:
        share_y_tl = check_bool_or_def(A1_opt[5], True)

    # Parse bonus options
    if "no_min_x_tl" in ut.A1_bonus:
        no_min_x_tl = True
    else:
        no_min_x_tl = False
    if "no_min_x_tl" in ut.A1_bonus:
        no_min_y_tl = True
    else:
        no_min_y_tl = False
    if "hide_shared_x_t" in ut.A1_bonus:
        hide_shared_x_t = True
    else:
        hide_shared_x_t = False
    if "hide_shared_y_t" in ut.A1_bonus:
        hide_shared_y_t = True
    else:
        hide_shared_y_t = False
    figsize = ut.extract_bonus("figsize", array_type=float)
    dpi = ut.extract_bonus("dpi", type=int, default=400)
    panel_labels = ut.extract_bonus("panel_labels")
    ix_match_y_lim = ut.extract_bonus("ix_match_y_lim", type=int)
    iy_match_x_lim = ut.extract_bonus("iy_match_x_lim", type=int)

    # Set up the figure
    if figsize is None:
        figsize = (9 + 7 * (num_x - 1), 9 + 7 * (num_y - 1))
    elif type(figsize) != tuple:
        figsize = tuple(figsize)
    fig = plt.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(nrows=num_y, ncols=num_x)
    A1_ax = []
    A2_ix_iy = []
    for i_y in range(num_y):
        for i_x in range(num_x):
            A1_ax.append(plt.subplot(gs[i_y, i_x]))
            A2_ix_iy.append([i_x, i_y])

    # ========
    # Plot each panel
    # ========
    for i_ax, (ax_i, A1_ix_iy, A1_arg) in enumerate(zip(A1_ax, A2_ix_iy, A2_arg)):
        plt.sca(ax_i)
        func_name = A1_arg[0]

        # Auto axis tweaks
        if func_name == "plot_snapshots":
            ax_i.set_aspect("equal")

        # Add any default "." arguments to fill up until A1_fig_ax
        while len(A1_arg) < len(inspect.getfullargspec(eval(func_name)).args):
            A1_arg.append(".")

        # Run the function to make the plot
        cmd = func_name + '("' + '", "'.join(A1_arg[1:]) + '", [fig, ax_i])'
        print("\nax_%d: %s" % (i_ax, cmd))
        exec(cmd)

        # Remove x label
        if share_x_l and A1_ix_iy[1] != num_y - 1:
            ax_i.set_xlabel("")
        # Remove y label
        if share_y_l and A1_ix_iy[0] != 0:
            ax_i.set_ylabel("")

        # Remove minimum x tick label
        if no_min_x_tl and A1_ix_iy[0] != 0 and A1_ix_iy[1] == num_y - 1:
            A1_label = ax_i.get_xticklabels()
            A1_label[0] = ""
            ax_i.set_xticklabels(A1_label)
        # Remove minimum y tick label
        if no_min_y_tl and A1_ix_iy[1] != num_y - 1:
            A1_label = ax_i.get_yticklabels()
            A1_label[0] = ""
            ax_i.set_yticklabels(A1_label)

        # Remove x tick labels (and ticks)
        if share_x_tl and A1_ix_iy[1] != num_y - 1:
            ax_i.set_xticklabels([])
            if ax_i.get_xscale() == "log":
                ax_i.set_xticklabels([], minor=True)
            if hide_shared_x_t:
                ax_i.tick_params("x", bottom=False, which="both")
        # Remove y tick labels (and ticks)
        if share_y_tl and A1_ix_iy[0] != 0:
            ax_i.set_yticklabels([])
            if ax_i.get_yscale() == "log":
                ax_i.set_yticklabels([], minor=True)
            if hide_shared_y_t:
                ax_i.tick_params("y", left=False, which="both")

        # Select axis limits to match
        if ix_match_y_lim is not None and A1_ix_iy[0] == ix_match_y_lim:
            y_lim_0, y_lim_1 = ax_i.get_ylim()
        if iy_match_x_lim is not None and A1_ix_iy[1] == iy_match_x_lim:
            x_lim_0, x_lim_1 = ax_i.get_xlim()

        # Shorten log major ticks
        if share_x_tl and A1_ix_iy[1] != num_y - 1 and ax_i.get_xscale() == "log":
            ax_i.xaxis.set_tick_params(width=1, length=6, which="major")
        if share_y_tl and A1_ix_iy[0] != 0 and ax_i.get_yscale() == "log":
            ax_i.yaxis.set_tick_params(width=1, length=6, which="major")

        # Panel labels
        if panel_labels is not None:
            # Parse
            v_loc = panel_labels[0]
            h_loc = panel_labels[1]
            label_type = panel_labels[2]

            # Location
            d = 0.03
            if v_loc == "u":
                y = 1 - d
                va = "top"
            elif v_loc == "l":
                y = d
                va = "bottom"
            else:
                raise Exception("Unexpected element in panel_labels ", panel_labels)
            if h_loc == "l":
                x = d
                ha = "left"
            elif h_loc == "c":
                x = 0.5
                ha = "left"
            elif h_loc == "r":
                x = 1 - d
                ha = "right"
            else:
                raise Exception("Unexpected element in panel_labels ", panel_labels)

            # Label
            if label_type == "a":
                A1_text = A1_enum_a[:num_tile]
            elif label_type == "A":
                A1_text = A1_enum_A[:num_tile]
            elif label_type == "1":
                A1_text = A1_enum_1[:num_tile]
            elif label_type == "i":
                A1_text = A1_enum_i[:num_tile]
            elif label_type == "I":
                A1_text = A1_enum_I[:num_tile]
            else:
                raise Exception("Unexpected element in panel_labels ", panel_labels)

            txt = ax_i.text(
                x,
                y,
                r"\textbf{%s}" % A1_text[i_ax],
                transform=ax_i.transAxes,
                fontsize=28 + 4 * (num_x - 2),
                color="k",
                ha=ha,
                va=va,
            )

            # Background outline
            txt.set_path_effects([path_effects.withStroke(linewidth=4, foreground="w")])

        print("ax_%d: Done" % (i_ax))

    # Set matching axis limits
    if ix_match_y_lim is not None or iy_match_x_lim is not None:
        for i_ax, (ax_i, A1_ix_iy, A1_arg) in enumerate(zip(A1_ax, A2_ix_iy, A2_arg)):
            if share_y_tl and ix_match_y_lim is not None:
                ax_i.set_ylim(y_lim_0, y_lim_1)
            if share_x_tl and iy_match_x_lim is not None:
                ax_i.set_xlim(x_lim_0, x_lim_1)

    # Adjust spacing
    if x_space is not None:
        plt.subplots_adjust(wspace=x_space)
    if y_space is not None:
        plt.subplots_adjust(hspace=y_space)

    # Save
    savefig(check_end(Fp_save, ".pdf"), dpi=dpi)
    print('\nSaved "%s" ' % Fp_save[-64:])
    plt.close()


def print_results_table(sim_set):
    """Print a latex table of results, with custom options for different cases.

    Redirect or tee to a txt file for easier copy pasting.

    Parameters
    ----------
    sim_set : SimSet
        An object with information about a set of simulations.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim_set, SimSet):
        sim_set = Di_sim_set[sim_set]
    sim_set.init()

    # Track simulation names and times to avoid duplicates
    A1_name_time_done = []
    if sim_set.category == "phodei":
        sim_set.make_unique()

        # SPH table list for reference in the REBOUND table
        if isinstance(sim_set.A1_sim[0], cl.RebSim):
            ob.set_phodei_tables.make_unique()
            A1_sim_in = np.array([sim.name for sim in ob.set_phodei_tables.A1_sim])

    # ========
    # Load simulations and extract data
    # ========
    i_sim = 0
    while i_sim < len(sim_set.A1_sim):
        sim = sim_set.A1_sim[i_sim]
        sim.init()

        # Snapshot time, automatically or from the set
        if isinstance(sim, cl.RebSim):
            rsim = sim
            time = -1
        elif sim.category == "phodei":
            time = sim.A1_time_snap[-1]
        else:
            time = int(sim_set.A1_time[i_sim])
        if isinstance(sim, cl.Simulation):
            time = sim.select_snapshot_times(sim_set.A1_time[i_sim])[0]

        # Avoid duplicate simulations
        if "%s %d" % (sim.name, time) in A1_name_time_done:
            i_sim += 1
            continue
        A1_name_time_done.append("%s %d" % (sim.name, time))

        # Size for the latex table(s), default to single column and no row splitting
        n_col = 1
        n_row = None

        # Load or calculate project-specific data
        if sim.category == "phodei":
            if isinstance(sim, cl.Simulation):
                if "A2000c30" in sim.name:
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
                    ) = um.accum_data_phodei_m_capt(sim)
                else:
                    (
                        m_bnd_in_Hill,
                        m_bnd_in_05Hill,
                        m_bnd_out_Hill,
                        m_unb,
                        m_bnd_in_Hill_fof_min,
                        num_fof_min,
                        num_fof_min_capt,
                    ) = um.accum_data_phodei_m_capt(sim)
            else:
                (
                    m_tot_in,
                    m_rm_close,
                    m_rm_far,
                    m_bump,
                    m_disr,
                    m_no_col,
                    m_surv,
                    m_col,
                ) = um.accum_data_phodei_reb_m(sim)

        # ========
        # Latex table entries
        # ========
        if i_sim == 0:
            A2_table = [[]]

        if sim.category == "phodei":
            if isinstance(sim, cl.Simulation):
                # Mass entry
                M_entry = "$10^{%.3g}$" % np.round(np.log10(sim.impact.M_2), 2)
                if "A2000c30" in sim.name:
                    M_entry += " d"

                # Spin entry
                if "spin=z" in sim.impact.A1_preset_2:
                    L_entry = "$z$ & $"
                elif "spin=-z" in sim.impact.A1_preset_2:
                    L_entry = "$z$ & $-"
                elif "spin=x" in sim.impact.A1_preset_2:
                    L_entry = "$x$ & $"
                elif "spin=y" in sim.impact.A1_preset_2:
                    L_entry = "$y$ & $"
                else:
                    L_entry = " & "
                if "_s030" in sim.name:
                    L_entry += "1$"
                elif "_s036" in sim.name:
                    L_entry += "3/4$"
                elif "_s047" in sim.name:
                    L_entry += "1/2$"
                elif "_s086" in sim.name:
                    L_entry += "1/4$"
                elif "_s170" in sim.name:
                    L_entry += "1/8$"

                # Store all entries
                entry = (
                    "    %d " % (i_sim + 1)
                    + "& $%.1f$ " % (sim.impact.r_p / Ma.R)
                    + "& $%.1f$ " % (sim.impact.v_inf / 1e3)
                    + "& %s " % M_entry
                    + "& %s " % sim.n_str
                    + "& %s " % L_entry
                    + "& %d " % num_fof_min
                    + "& $%.1f$ " % (m_unb / sim.impact.M_2 * 100)
                    + "& $%.1f$ "
                    % ((m_bnd_in_Hill + m_bnd_out_Hill) / sim.impact.M_2 * 100)
                    + "& $%.1f$ " % (m_bnd_in_Hill / sim.impact.M_2 * 100)
                    + "& $%.1f$ " % (m_bnd_in_Hill_fof_min / sim.impact.M_2 * 100)
                    + "& %d " % num_fof_min_capt
                )
            else:
                # Store all entries
                entry = (
                    "    %d " % (np.where(A1_sim_in == sim.sim_in.name)[0][0] + 1)
                    + "& $%d$ " % (sim.sim_in_angle_2)
                    + "& $%d$ " % (-sim.sim_in_angle)
                    + "& $%d$ " % (-sim.oblate_obliquity)
                    + "& %s "
                    % (
                        "$%d$ " % (sim.A1_picle[1].nu * rad_to_deg)
                        if "nu" in sim.name
                        else ""
                    )
                    + "& $%.2f$ " % (m_rm_close / 1e19)
                    + "& $%.2f$ " % (m_rm_far / 1e19)
                    + "& $%.2f$ " % (m_bump / 1e19)
                    + "& $%.2f$ " % (m_surv / 1e19)
                    + "& $%.2f$ " % (m_col / 1e19)
                )

            A2_table[0].append(entry)
        else:
            raise Exception(
                "Latex table not yet implemented for set, category ",
                sim_set.name,
                ",",
                sim.category,
            )

        # Next simulation
        i_sim += 1

    # ========
    # Print
    # ========
    print("\n\n# ========\n")
    # Split into separate subtables if too long
    if n_row is not None:
        while n_row * n_col < len(A2_table[-1]):
            # Add the next subtable
            A2_table.append(A2_table[-1][n_row * n_col :])
            # Remove from the previous subtable
            A2_table[-2] = A2_table[-2][: n_row * n_col]

    # Print each subtable
    for i_tab in range(len(A2_table)):
        # Split the subtable into columns
        A2_subtable = np.array_split(A2_table[i_tab], n_col)

        # Print each row of the subtable
        for row in range(len(A2_subtable[0])):
            try:
                if row != 0:
                    # Indent all but the first row for easy collapsing in the latex editor
                    print("  ", end="")

                for col in range(n_col - 1):
                    print(A2_subtable[col][row], end=" & \n& ")
                print(A2_subtable[-1][row], end=" \\\\ \n")
            except IndexError:
                print("\\\\")

        # Separate from the next subtable
        print("\n")


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR custom_plots.py  ====\n")
