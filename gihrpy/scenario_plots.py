"""
Functions for GIHR plots of accumulated scenario results, etc.
"""

from jkeger import *
import gihrpy.utilities as ut
from gihrpy.classes_sim import Simulation, SimSet
from gihrpy.objects import Di_sim, Di_sim_set


def plot_sim_param_subsets(
    sim_set,
    param_y,
    param_x,
    param_c=None,
    param_ls=None,
    time_all=None,
    A1_fig_ax=None,
):
    """Plot a result parameter from each simulation in a set.

    Optionally divided into subsets by discrete parameters, distinguished by
    colour and/or by marker or linestyle.

    Plot as a function of the chosen x-axis parameter, either as points or lines.

    e.g. captured mass as a function of encounter periapsis, in subsets by spin.

    Parameters
    ----------
    sim_set : SimSet
        An object with information about a set of simulations.

    param_y, param_x : str
        The parameters to plot on the y and x axes. See Simulation.param_value.

    param_c, param_ls : str (opt.)
        Split the data into subsets by these parameters, plotted by colour
        and/or by linestyle. See Simulation.param_value.

    time_all : int
        The time of the snapshot to plot (s), if not already set by sim_set.
        0       The initial snapshot.
        -1      The final snapshot.

    A1_fig_ax : [fig, ax]
        If provided, then plot on this existing figure and axes instead of
        making new ones.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim_set, SimSet):
        sim_set = Di_sim_set[sim_set]
    sim_set.init()
    po = sim_set.po

    param_c = check_none_or_def(param_c, None)
    param_ls = check_none_or_def(param_ls, None)
    time_all = check_none_or_def(time_all, None)

    # Figure
    if A1_fig_ax is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        fig, ax = A1_fig_ax
        plt.figure(fig.number)

    # Subset parameter values : colours/linestyles must be set in advance
    if param_c is not None:
        Di_p_colour = sim_set.Di_param_Di_p_colour[param_c]
    if param_ls is not None:
        Di_p_linestyle = sim_set.Di_param_Di_p_linestyle[param_ls]

    # Track used param values
    A1_p_c_all = []
    A1_p_ls_all = []

    # Initialise arrays for each unique p_c, p_ls combination
    def id_c_ls(p_c, p_ls):
        """Simple string to identify a subset parameter combination."""
        try:
            return "%g %g" % (p_c, p_ls)
        except TypeError:
            try:
                return "%s %g" % (p_c, p_ls)
            except TypeError:
                try:
                    return "%g %s" % (p_c, p_ls)
                except TypeError:
                    return "%s %s" % (p_c, p_ls)

    Di_id_A1_p_y = {
        id_c_ls(p_c, p_ls): []
        for p_ls in Di_p_linestyle.keys()
        for p_c in Di_p_colour.keys()
    }
    Di_id_A1_p_x = deepcopy(Di_id_A1_p_y)
    A1_p_x_all = []
    A1_p_y_all = []
    # Custom arrays for other values
    if param_y == "m_capt" and "A2000c30" in sim_set.name:
        Di_id_A1_m_c_capt = deepcopy(Di_id_A1_p_y)

    # Ensure each simulation is only loaded once
    sim_set.ensure_unique()

    # ========
    # Load each simulation and store data
    # ========
    for i_sim, sim in enumerate(sim_set.A1_sim):
        sim.init()
        po = sim.po

        # Parameter values (plotting units)
        p_x = sim.param_value(param_x)
        p_y = sim.param_value(param_y)
        if param_c is not None:
            p_c = sim.param_value(param_c)
        else:
            p_c = 0
        if param_ls is not None:
            p_ls = sim.param_value(param_ls)
        else:
            p_ls = 0

        # Store values
        Di_id_A1_p_x[id_c_ls(p_c, p_ls)].append(p_x)
        Di_id_A1_p_y[id_c_ls(p_c, p_ls)].append(p_y)
        A1_p_x_all.append(p_x)
        A1_p_y_all.append(p_y)
        A1_p_c_all.append(p_c)
        A1_p_ls_all.append(p_ls)

        # Store other values in custom arrays for later plotting or fitting
        if param_y == "m_capt" and "A2000c30" in sim_set.name:
            # Captured mass fraction of core material
            if "A2000c30" in sim.name:
                Di_id_A1_m_c_capt[id_c_ls(p_c, p_ls)].append(
                    sim.param_value("m_c_capt")
                )

    # Reference simulation for label utilities etc
    sim = sim_set.A1_sim[0]
    po = sim_set.po

    print("")
    A1_p_x_all = np.array(A1_p_x_all)
    A1_p_y_all = np.array(A1_p_y_all)
    A1_p_c_all = np.array(A1_p_c_all)
    A1_p_ls_all = np.array(A1_p_ls_all)

    # ========
    # Plot the lines
    # ========
    for i_ls, p_ls in enumerate(Di_p_linestyle.keys()):
        for i_c, p_c in enumerate(Di_p_colour.keys()):
            # Extract the data in this subset
            id_i = id_c_ls(p_c, p_ls)
            A1_p_x = np.array(Di_id_A1_p_x[id_i])
            A1_p_y = np.array(Di_id_A1_p_y[id_i])

            # Skip empty subsets
            if len(A1_p_x) == 0:
                continue

            # Ensure the right order
            A1_sort_x = np.argsort(A1_p_x)

            # Plot
            ax.plot(
                A1_p_x[A1_sort_x],
                A1_p_y[A1_sort_x],
                c=Di_p_colour[p_c],
                ls=Di_p_linestyle[p_ls],
                alpha=1,
            )

            # Plot other values
            if param_y == "m_capt" and "A2000c30" in sim_set.name:
                ax.plot(
                    A1_p_x[A1_sort_x],
                    np.array(Di_id_A1_m_c_capt[id_i])[A1_sort_x],
                    c=Di_p_colour[p_c],
                    ls=ls_dash,
                    alpha=alpha,
                )

    # Remove unused elements from the dictionaries for the legends
    for p_c in [key for key in Di_p_colour.keys()]:
        if p_c not in A1_p_c_all:
            del Di_p_colour[p_c]
    for p_ls in [key for key in Di_p_linestyle.keys()]:
        if p_ls not in A1_p_ls_all:
            del Di_p_linestyle[p_ls]

    # Default legend parameters
    loc_p_c = "lower right"
    loc_p_ls = "upper left"
    c_leg = "0.2"
    fs_leg = 22
    ncol_p_c = 1
    ncol_p_ls = 1
    num_p_c = len(np.unique(A1_p_c_all))
    num_p_ls = len(np.unique(A1_p_ls_all))

    # Legend custom locations or remove
    if param_y == "m_capt":
        loc_p_c = "lower left"
    if sim_set.name == "set_Ma_xp_A2000_s___z_n65_r___v__":
        loc_p_c = "lower right"
        loc_p_ls = "upper right"
    if param_c == "_r":
        # No legend for reoriented repeats
        num_p_c = 0
    if param_ls == "_r":
        # No legend for reoriented repeats
        num_p_ls = 0

    # Colour legend
    if num_p_c > 1:
        # Legend entries
        A1_leg_p_c = []
        for p_c, c in Di_p_colour.items():
            A1_leg_p_c.append(
                ax.plot([], [], color=c, label=sim.param_value_str(param_c, value=p_c))
            )

        # Flatten list if needed
        if isinstance(A1_leg_p_c[0], list):
            A1_leg_p_c = [item for sublist in A1_leg_p_c for item in sublist]

        # Add legend
        leg_p_c = ax.legend(
            handles=A1_leg_p_c,
            prop={"size": fs_leg},
            loc=loc_p_c,
            title=sim.param_label_unit(param_c),
            ncol=ncol_p_c,
        )
        plt.setp(leg_p_c.get_title(), fontsize=fs_leg + 2)
        ax.add_artist(leg_p_c)

    # Linestyle legend
    if num_p_ls > 1:
        # Legend entries
        A1_leg_p_ls = []
        for p_ls, ls in Di_p_linestyle.items():
            # Skip unwanted entries
            if sim_set.name == "set_Ma_xp_A2000_s___z_n65_r___v__" and p_ls == -1:
                continue

            A1_leg_p_ls.append(
                ax.plot(
                    [],
                    [],
                    c="k",
                    ls=ls,
                    label=sim.param_value_str(param_ls, value=p_ls),
                )
            )

        # Flatten list if needed
        if isinstance(A1_leg_p_ls[0], list):
            A1_leg_p_ls = [item for sublist in A1_leg_p_ls for item in sublist]

        # Add legend
        leg_p_ls = ax.legend(
            handles=A1_leg_p_ls,
            prop={"size": fs_leg},
            loc=loc_p_ls,
            title=sim.param_label_unit(param_ls),
            ncol=ncol_p_ls,
        )
        plt.setp(leg_p_ls.get_title(), fontsize=fs_leg + 2)
        ax.add_artist(leg_p_ls)

    # ========
    # Axes etc
    # ========
    ax.set_xlabel(sim.param_label_unit(param_x))
    ax.set_ylabel(sim.param_label_unit(param_y))

    if param_x in po.Di_param_log.keys() and po.Di_param_log[param_x]:
        ax.set_xscale("log")
    if param_y in po.Di_param_log.keys() and po.Di_param_log[param_y]:
        ax.set_yscale("log")

    if param_x in po.Di_param_lim.keys():
        x_lim_0, x_lim_1 = po.Di_param_lim[param_x]
        ax.set_xlim(x_lim_0, x_lim_1)
    if param_y in po.Di_param_lim.keys():
        y_lim_0, y_lim_1 = po.Di_param_lim[param_y]
        ax.set_ylim(y_lim_0, y_lim_1)

    # Misc extras
    if po.func_param_subsets_extra is not None:
        po.func_param_subsets_extra(
            sim_set,
            param_x,
            param_y,
            param_c,
            param_ls,
            ax,
        )

    auto_nice_plot(ax, fontsize=po.fontsize)

    # Save and close, unless adding to an existing figure
    if A1_fig_ax is None:
        Fp_save = "%s/%s_subs_%s_%s" % (
            po.dir_save,
            sim_set.name,
            param_y,
            param_x,
        )
        if num_p_c > 1:
            Fp_save += "_%s" % param_c
        if num_p_ls > 1:
            Fp_save += "_%s" % param_ls
        if time_all is not None:
            Fp_save += "_%06d" % sim.select_snapshot_times(time_all)[0]
        if po.Fp_suffix is not None:
            Fp_save += "%s" % po.Fp_suffix

        # Save
        if ut.do_paper:
            savefig(Fp_save + ".pdf", dpi=400)
        else:
            savefig(Fp_save + ".png", dpi=po.dpi)
        print('Saved "%s" ' % Fp_save[-64:])

        plt.close()


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR scenario_plots.py  ====\n")
