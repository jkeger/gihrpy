"""
Functions for GIHR initial conditions, etc.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.utilities_reb as ur
from gihrpy.classes import PlotOptions, InitProf, InitCond
from gihrpy.classes_sim import Simulation
from gihrpy.classes_reb import RebSim
from gihrpy.objects import Di_init_prof, Di_init_cond, Di_sim, Di_reb_sim


def _plot_profiles(
    A1_r,
    A1_rho,
    A1_u,
    A1_P,
    A1_T,
    A1_mat_id,
    A1_r_picle=None,
    A1_rho_picle=None,
    A1_u_picle=None,
    A1_P_picle=None,
    A1_T_picle=None,
    A1_mat_id_picle=None,
    po=PlotOptions(),
    P_unit=None,
    P_label=None,
    c=None,
    ls="-",
    alpha=0.8,
    A1_fig_ax=None,
):
    """Plot a planetary profile, and particles if provided.

    Parameters
    ----------
    A1_r : [float]
        The profile radii (m).

    A1_rho : [float]
        The profile densities (kg m^-3).

    A1_u : [float]
        The profile specific internal energies (J kg^-1).

    A1_P : [float]
        The profile pressures (Pa).

    A1_T : [float]
        The profile temperatures (K).

    A1_mat_id : [int]
        The profile material IDs.

    A1_*_picle : [float] or [int] (opt.)
        The same properties as the profiles above but for a set of particles.

    po : PlotOptions
        An object with plotting option settings.

    P_unit, P_label : float, str
        The units and units label for the pressures, e.g. 1e9 and "GPa".

    c, ls, alpha : str
        A colour, linestyle, and transparency to plot the profile lines.
        Defaults to colour each material by their colour in po.Di_mat_colour.

    A1_fig_ax : [fig, ax1, ax2, ax3, ax4]
        If provided, then plot on this existing figure and axes instead of
        making new ones.

    Returns
    -------
    A1_fig_ax : [fig, ax1, ax2, ax3, ax4]
        The figure and axes.
    """
    # Particles or not
    if A1_r_picle is not None:
        do_picle = True
    else:
        do_picle = False

    # Plotting units
    A1_r /= po.r_unit
    if do_picle:
        A1_r_picle /= po.r_unit
    if P_unit is None:
        if np.amax(A1_P) > 1e8:
            P_unit = 1e9
            P_label = "GPa"
        else:
            P_unit = 1e6
            P_label = "MPa"

    # Index of the start of each layer
    A1_idx_layer = np.concatenate(
        ([0], np.where(np.diff(A1_mat_id) != 0)[0] + 1, [len(A1_r)])
    )

    # Colour for each layer
    if c is None:
        A1_c_layer = [
            po.Di_id_colour[mat_id] for mat_id in A1_mat_id[A1_idx_layer[:-1]]
        ]
        if do_picle:
            # Darken profile colours
            A1_c_layer = [tweak_luminosity(c, 0.4) for c in A1_c_layer]

            # Particle colours
            A1_colour_picle = np.empty(len(A1_r_picle), dtype=object)
            for mat_id, colour in po.Di_id_colour.items():
                A1_sel = np.where(A1_mat_id_picle == mat_id)[0]
                A1_colour_picle[A1_sel] = colour
    else:
        A1_c_layer = [c] * len(A1_idx_layer)
        if do_picle:
            A1_colour_picle = c

    # Figure
    if A1_fig_ax is None:
        fig = plt.figure(figsize=(9, 9))
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=2)
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[0, 1])
        ax3 = plt.subplot(gs[1, 0])
        ax4 = plt.subplot(gs[1, 1])
        gs.update(wspace=0.02, hspace=0.02)
    else:
        fig, ax1, ax2, ax3, ax4 = A1_fig_ax
        plt.figure(fig.number)

    # Rasterise the particles (anything with zorder < 1)
    for ax in fig.axes:
        ax.set_rasterization_zorder(1)

    # Profiles
    if do_picle:
        alpha *= 0.7
    for ax, A1_param in zip(
        [ax1, ax2, ax3, ax4], [A1_rho, A1_u / 1e6, A1_P / P_unit, A1_T]
    ):
        for i_layer in range(len(A1_idx_layer) - 1):
            ax.plot(
                A1_r[A1_idx_layer[i_layer] : A1_idx_layer[i_layer + 1]],
                A1_param[A1_idx_layer[i_layer] : A1_idx_layer[i_layer + 1]],
                c=A1_c_layer[i_layer],
                ls=ls,
                alpha=alpha,
                zorder=2,
            )

    # Particles
    if do_picle:
        size = 2**2
        alpha *= 0.4
        for ax, A1_param in zip(
            [ax1, ax2, ax3, ax4],
            [A1_rho_picle, A1_u_picle / 1e6, A1_P_picle / P_unit, A1_T_picle],
        ):
            ax.scatter(
                A1_r_picle,
                A1_param,
                c=A1_colour_picle,
                edgecolor="none",
                s=size,
                marker=".",
                alpha=alpha,
                zorder=0,
            )

    # Axes etc.
    ax1.set_ylabel(r"Density (kg m$^{-1}$)")
    ax2.set_ylabel(r"Sp. Int. Energy (MJ kg$^{-1}$)")
    ax3.set_ylabel(r"Pressure (%s)" % P_label)
    ax4.set_ylabel(r"Temperature (K)")
    for ax in [ax1, ax2]:
        ax.set_xticklabels([])
        ax.set_xticklabels([], minor=True)
    for ax in [ax3, ax4]:
        ax.set_xlabel(r"Radius (%s)" % po.R_label)
    for ax in [ax2, ax4]:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim(0, None)

    auto_nice_plot(do_far_ticks=False)

    return [fig, ax1, ax2, ax3, ax4]


def gen_init_prof(ip, check_min_period=None, verbosity=None):
    """Generate a planetary profile for initial conditions, using WoMa.

    Parameters
    ----------
    ip : InitProf
        The planet profile object.

    check_min_period : bool
        For spinning planets, checks if the period provided is less than the
        minimum physically allowed. Slow -- set True only if required for
        extremely high spin.

    verbosity : int
        Override for ip.verbosity.
    """
    # Set the object if the name was provided instead
    if not isinstance(ip, InitProf):
        ip = Di_init_prof[ip]
    ip.init()

    verbosity = check_none_or_def(verbosity, None)
    if verbosity is not None:
        ip.verbosity = int(verbosity)
    check_min_period = check_bool_or_def(check_min_period, False)
    woma.load_eos_tables(ip.A1_mat_layer)

    os.makedirs(os.path.dirname(ip.Fp_save), exist_ok=True)

    # ========
    # Make and save the profiles
    # ========
    if ip.period is None:
        # Call the appropriate function for the given number of layers and radius/mass inputs
        if ip.num_layer == 1:
            print("  gen_prof_L1_find_R_given_M()")
            ip.gen_prof_L1_find_R_given_M(
                R_max=ip.R_max,
                tol=ip.tol,
                tol_M_tweak=ip.tol_M_tweak,
                num_attempt=ip.num_attempt,
                verbosity=ip.verbosity,
            )
        elif ip.num_layer == 2:
            if ip.R_min is not None:
                print("  gen_prof_L2_find_R_R1_given_M1_M2()")
                ip.gen_prof_L2_find_R_R1_given_M1_M2(
                    R_min=ip.R_min,
                    R_max=ip.R_max,
                    tol=ip.tol,
                    tol_M_tweak=ip.tol_M_tweak,
                    num_attempt=ip.num_attempt,
                    verbosity=ip.verbosity,
                )
            else:
                print("  gen_prof_L2_find_R1_given_M_R()")
                ip.gen_prof_L2_find_R1_given_M_R(
                    tol=ip.tol,
                    tol_M_tweak=ip.tol_M_tweak,
                    num_attempt=ip.num_attempt,
                    verbosity=ip.verbosity,
                )

        ip.save(ip.Fp_save)
    else:
        if check_min_period:
            ip.verbosity += 1

        # Load the base spherical planet
        planet = woma.Planet(load_file=ip.planet.Fp_save)

        # Make the profiles!
        spin_planet = woma.SpinPlanet(
            planet=planet,
            period=ip.period,
            num_prof=ip.num_prof_spin,
            f_iter=ip.f_iter,
            check_min_period=check_min_period,
            tol_density_profile=ip.tol_density_profile,
            tol_layer_masses=ip.tol_layer_masses,
            num_attempt_1=ip.num_attempt_1,
            num_attempt_2=ip.num_attempt_2,
            verbosity=ip.verbosity,
        )

        spin_planet.save(ip.Fp_save)

    # Units
    m_print = ip.po.Di_param_unit_label["m"].unit_print
    r_print = ip.po.Di_param_unit_label["r"].unit_print

    # Print outputs for future reference
    print("\nAttributes to copy-paste into %s's declaration for reference:  " % ip.name)
    space = 15
    if ip.period is None:
        print(
            "        %s = np.array(%s)%s,"
            % (
                add_whitespace("A1_M_layer", space),
                format_array_string(ip.A1_M_layer / ip.po.M_unit, "%.5g"),
                m_print,
            )
        )
        print(
            "        %s = np.array(%s)%s,"
            % (
                add_whitespace("A1_R_layer", space),
                format_array_string(ip.A1_R_layer / ip.po.r_unit, "%.5g"),
                r_print,
            )
        )
        print(
            "        %s = %s,"
            % (
                add_whitespace("A1_idx_layer", space),
                format_array_string(ip.A1_idx_layer, "%d"),
            )
        )
        print(
            "        %s = %.5g%s,"
            % (add_whitespace("M", space), ip.M / ip.po.M_unit, m_print)
        )
        print("        %s = %.5g," % (add_whitespace("P_s", space), ip.P_s))
        print("        %s = %.5g," % (add_whitespace("T_s", space), ip.T_s))
        print("        %s = %.5g," % (add_whitespace("rho_s", space), ip.rho_s))
        if ip.num_layer > 2:
            print("        %s = %.5g," % (add_whitespace("P_2", space), ip.P_2))
            print("        %s = %.5g," % (add_whitespace("T_2", space), ip.T_2))
            print("        %s = %.5g," % (add_whitespace("rho_2", space), ip.rho_2))
        if ip.num_layer > 1:
            print("        %s = %.5g," % (add_whitespace("P_1", space), ip.P_1))
            print("        %s = %.5g," % (add_whitespace("T_1", space), ip.T_1))
            print("        %s = %.5g," % (add_whitespace("rho_1", space), ip.rho_1))
        print("        %s = %.5g," % (add_whitespace("P_0", space), ip.P_0))
        print("        %s = %.5g," % (add_whitespace("T_0", space), ip.T_0))
        print("        %s = %.5g," % (add_whitespace("rho_0", space), ip.rho_0))
        print("        %s = %.5g," % (add_whitespace("I_MR2", space), ip.I_MR2))
    else:
        print(
            "        %s = %.5g," % (add_whitespace("period", space), spin_planet.period)
        )
        print(
            "        %s = np.array(%s)%s,"
            % (
                add_whitespace("A1_M_layer", space),
                format_array_string(spin_planet.A1_M_layer / ip.po.M_unit, "%.5g"),
                m_print,
            )
        )
        print(
            "        %s = np.array(%s)%s,"
            % (
                add_whitespace("A1_R_layer", space),
                format_array_string(spin_planet.A1_R_layer / ip.po.r_unit, "%.5g"),
                r_print,
            )
        )
        print(
            "        %s = np.array(%s)%s,"
            % (
                add_whitespace("A1_Z_layer", space),
                format_array_string(spin_planet.A1_Z_layer / ip.po.r_unit, "%.5g"),
                r_print,
            )
        )
        print(
            "        %s = %s,"
            % (
                add_whitespace("A1_idx_layer", space),
                format_array_string(spin_planet.A1_idx_layer_eq, "%d"),
            )
        )
        print(
            "        %s = %.5g%s,"
            % (add_whitespace("M", space), spin_planet.M / ip.po.M_unit, m_print)
        )
        print("        %s = %.5g," % (add_whitespace("P_s", space), spin_planet.P_s))
        print("        %s = %.5g," % (add_whitespace("T_s", space), spin_planet.T_s))
        print(
            "        %s = %.5g," % (add_whitespace("rho_s", space), spin_planet.rho_s)
        )
        if spin_planet.num_layer > 2:
            print(
                "        %s = %.5g," % (add_whitespace("P_2", space), spin_planet.P_2)
            )
            print(
                "        %s = %.5g," % (add_whitespace("T_2", space), spin_planet.T_2)
            )
            print(
                "        %s = %.5g,"
                % (add_whitespace("rho_2", space), spin_planet.rho_2)
            )
        if spin_planet.num_layer > 1:
            print(
                "        %s = %.5g," % (add_whitespace("P_1", space), spin_planet.P_1)
            )
            print(
                "        %s = %.5g," % (add_whitespace("T_1", space), spin_planet.T_1)
            )
            print(
                "        %s = %.5g,"
                % (add_whitespace("rho_1", space), spin_planet.rho_1)
            )
        print("        %s = %.5g," % (add_whitespace("P_0", space), spin_planet.P_0))
        print("        %s = %.5g," % (add_whitespace("T_0", space), spin_planet.T_0))
        print(
            "        %s = %.5g," % (add_whitespace("rho_0", space), spin_planet.rho_0)
        )
        print(
            "        %s = %.5g," % (add_whitespace("I_MR2", space), spin_planet.I_MR2)
        )
        print(
            "        %s = %.5g,  # %.3g L_EM"
            % (add_whitespace("L", space), spin_planet.L, spin_planet.L / L_EM)
        )


def plot_init_prof(A1_ip, A1_colour=None, A1_label=None, Fp_save=None):
    """Plot one or a set of initial profiles with _plot_profiles().

    Parameters
    ----------
    A1_ip : InitProf or [InitProf]
        The name or a list of names of one or more profile objects to plot.
        Enter as a list including the brackets, commas, and no spaces, e.g.
        [ip_1,ip_2,ip_3]. See jkeger.array_from_string().

    A1_colour : [int] or "rbow"
        A list of colours for each profile set by an index for jkeger.A1_c, e.g.
        [0,1,2], or "rbow" for automatic rainbow colours in the order of A1_ip.

    A1_label : [str] or "name"
        A list of labels for each profile, or "name" to use the profile names.

    Fp_save : str (opt.)
        The filename for saving the figure. Defaults to the combined names.
    """
    A1_colour = check_none_or_def(A1_colour, None)
    A1_label = check_none_or_def(A1_label, None)

    # Set the objects from the names
    try:
        A1_ip = [Di_init_prof[A1_ip]]
    except KeyError:
        A1_ip = [Di_init_prof[ip] for ip in array_from_string(A1_ip, str)]
    num_ip = len(A1_ip)
    for ip in A1_ip:
        ip.init()

    # Filename
    if Fp_save is None:
        Fp_save = "%s" % A1_ip[0].name
        if num_ip > 1:
            Fp_save = "set_%s" % Fp_save
            for ip in A1_ip[1:]:
                Fp_save += "_%s" % ip.name
    if ip.po.Fp_suffix is not None:
        Fp_save += "_%s" % ip.po.Fp_suffix
    Fp_save = "%s/init_cond/%s" % (ip.po.dir_save, Fp_save)
    # Skip if already done
    if ut.do_skip_exist:
        if os.path.isfile(Fp_save):
            return

    # Colours
    if A1_colour is None:
        if num_ip == 1:
            A1_colour = [None]
        else:
            A1_colour = A1_c[:num_ip]
    elif A1_colour == "rbow":
        A1_colour = cmap_rbow(np.linspace(0, 1, num_ip))
    else:
        A1_colour = np.array(A1_c)[array_from_string(A1_colour, int)]

    # Labels
    if A1_label is None:
        A1_label = [None] * num_ip
    elif A1_label == "name":
        A1_label = [no_latex_string(ip.name) for ip in A1_ip]
    else:
        A1_label = array_from_string(A1_label, str)

    # Load and plot each profile
    A1_fig_ax = None
    for i_ip, (ip, c, l) in enumerate(zip(A1_ip, A1_colour, A1_label)):
        # Load the profiles
        ip.load_profiles()

        if i_ip == 0:
            if np.amax(ip.A1_P) > 1e8:
                P_unit = 1e9
                P_label = "GPa"
            else:
                P_unit = 1e6
                P_label = "MPa"

        # Plot!
        if ip.period is None:
            A1_fig_ax = _plot_profiles(
                A1_r=ip.A1_r,
                A1_rho=ip.A1_rho,
                A1_u=ip.A1_u,
                A1_P=ip.A1_P,
                A1_T=ip.A1_T,
                A1_mat_id=ip.A1_mat_id,
                po=ip.po,
                P_unit=P_unit,
                P_label=P_label,
                c=c,
                A1_fig_ax=A1_fig_ax,
            )
        else:
            # Both equatorial and polar profiles
            A1_fig_ax = _plot_profiles(
                A1_r=ip.A1_R,
                A1_rho=ip.A1_rho,
                A1_u=ip.A1_u,
                A1_P=ip.A1_P,
                A1_T=ip.A1_T,
                A1_mat_id=ip.A1_mat_id,
                po=ip.po,
                P_unit=P_unit,
                P_label=P_label,
                c=c,
                A1_fig_ax=A1_fig_ax,
            )
            A1_fig_ax = _plot_profiles(
                A1_r=ip.A1_Z,
                A1_rho=ip.A1_rho,
                A1_u=ip.A1_u,
                A1_P=ip.A1_P,
                A1_T=ip.A1_T,
                A1_mat_id=ip.A1_mat_id,
                po=ip.po,
                P_unit=P_unit,
                P_label=P_label,
                c=c,
                A1_fig_ax=A1_fig_ax,
            )

        # Values for axis limits
        if i_ip == 0:
            if ip.period is None:
                r_min = np.amin(ip.A1_r)
                r_max = np.amax(ip.A1_r)
            else:
                r_min = np.amin(ip.A1_Z)
                r_max = np.amax(ip.A1_R)
            rho_min = np.amin(ip.A1_rho)
            rho_max = np.amax(ip.A1_rho)
            u_min = np.amin(ip.A1_u)
            u_max = np.amax(ip.A1_u)
            P_min = np.amin(ip.A1_P)
            P_max = np.amax(ip.A1_P)
            T_min = np.amin(ip.A1_T)
            T_max = np.amax(ip.A1_T)
        else:
            if ip.period is None:
                r_min = min(r_min, np.amin(ip.A1_r))
                r_max = max(r_max, np.amax(ip.A1_r))
            else:
                r_min = min(r_min, np.amin(ip.A1_Z))
                r_max = max(r_max, np.amax(ip.A1_R))
            rho_min = min(rho_min, np.amin(ip.A1_rho))
            rho_max = max(rho_max, np.amax(ip.A1_rho))
            u_min = min(u_min, np.amin(ip.A1_u))
            u_max = max(u_max, np.amax(ip.A1_u))
            P_min = min(P_min, np.amin(ip.A1_P))
            P_max = max(P_max, np.amax(ip.A1_P))
            T_min = min(T_min, np.amin(ip.A1_T))
            T_max = max(T_max, np.amax(ip.A1_T))

        # For the legend
        A1_fig_ax[-1].plot([], [], c=c, label=l)
    print("")

    # Axes etc.
    set_auto_limits(A1_fig_ax[1], [r_min, r_max], [rho_min, rho_max], min_x=0)
    set_auto_limits(A1_fig_ax[2], [r_min, r_max], [u_min / 1e6, u_max / 1e6], min_x=0)
    set_auto_limits(
        A1_fig_ax[3], [r_min, r_max], [P_min / P_unit, P_max / P_unit], min_x=0
    )
    set_auto_limits(A1_fig_ax[4], [r_min, r_max], [T_min, T_max], min_x=0, min_y=0)
    if A1_label[0] is not None:
        A1_fig_ax[-1].legend()

    tight_layout()

    # Save
    if ut.do_paper:
        savefig(Fp_save + ".pdf", dpi=400, bbox_inches="tight")
    else:
        savefig(Fp_save + ".png", dpi=ip.po.dpi, bbox_inches="tight")

    print('Saved "%s" ' % Fp_save[-64:])
    plt.close()


def gen_init_cond(ic):
    """Generate particle initial conditions ready for a settling simulation.

    Parameters
    ----------
    ic : InitCond
        The initial conditions object.
    """
    # Set the object if the name was provided instead
    if not isinstance(ic, InitCond):
        ic = Di_init_cond[ic]
    ic.init()

    # Load the profiles
    if ic.ip.period is None:
        ic.ip.load_profiles()
    else:
        # Needed to set the right class type for woma checks
        planet = ic.ip.planet
        ic.ip = woma.SpinPlanet(load_file=ic.ip.Fp_save)

    # Place the particles
    particles = woma.ParticlePlanet(ic.ip, ic.num_picle_des, ic.N_ngb, verbosity=1)
    print("")

    # Total mass
    M = sum(particles.A1_m)

    # Approx minimum separation (central tetrahedron)
    min_sep = 1.63 * np.amin(root_sum_sq(particles.A2_pos.T))

    # Entropies if available for this EoS
    if any(["Til_" in mat for mat in ic.ip.A1_mat_layer]):
        do_entropies = False
    else:
        do_entropies = True

    # Save the particles
    particles.save(
        ic.Fp_save,
        boxsize=ic.boxsize,
        file_to_SI=ic.file_to_SI,
        do_entropies=do_entropies,
        verbosity=1,
    )

    # Print info for easy copy-pasting into the InitCond object
    M_unit = ic.po.Di_param_unit_label["m"].unit
    M_code = ic.po.Di_param_unit_label["m"].unit_code
    r_unit = ic.po.Di_param_unit_label["r"].unit
    R_code = ic.po.Di_param_unit_label["r"].unit_code
    print("\nAttributes to copy-paste into %s's declaration:  " % ic.name)
    print("        num_picle       = %d, " % particles.N_particles)
    print("        M               = %.6g%s," % (M / M_unit, M_code))
    print("        R_s             = %.6g%s," % (ic.R_s / r_unit, R_code))
    if ic.ip.period is not None:
        print("        Z_s             = %.6g%s," % (ic.Z_s / r_unit, R_code))
        print("        R_s_sph         = %.6g%s," % (ic.R_s_sph / r_unit, R_code))
    print("        min_sep         = %.3g * 1e3," % (min_sep / 1e3))


def gen_impact_init_cond(sim):
    """Generate initial conditions for an impact or similar simulation.

    Parameters
    ----------
    sim : Simulation
        An object with all the information about the simulation.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim, Simulation):
        sim = Di_sim[sim]
    sim.init()

    im = sim.impact
    im.prep_scenario()

    # Particle data
    Fp_init_cond_1 = check_end(im.Fp_init_cond_1, ".hdf5")
    Fp_init_cond_2 = check_end(im.Fp_init_cond_2, ".hdf5")
    print('Combining (1) "%s" ' % Fp_init_cond_1[-56:])
    print('      and (2) "%s" ' % Fp_init_cond_2[-56:])

    # Extract a fof group instead of all particles
    do_extract_fof_1 = False
    do_extract_fof_2 = False
    # Change the file path to the converted snapshot that was used as the fof input
    if im.A1_preset_1 is not None and im.A1_preset_1[0][:4] == "fof=":
        Fp_init_cond_1 = Fp_init_cond_1[:-5] + "_init.hdf5"
        do_extract_fof_1 = True
    if im.A1_preset_2 is not None and im.A1_preset_2[0][:4] == "fof=":
        Fp_init_cond_2 = Fp_init_cond_2[:-5] + "_init.hdf5"
        do_extract_fof_2 = True

    # Load particle data
    (
        A2_pos_1,
        A2_vel_1,
        A1_m_1,
        A1_h_1,
        A1_rho_1,
        A1_P_1,
        A1_u_1,
        A1_mat_id_1,
    ) = ut.load_particle_data_Fp(
        Fp_init_cond_1, ["pos", "vel", "m", "h", "rho", "P", "u", "mat_id"]
    )
    (
        A2_pos_2,
        A2_vel_2,
        A1_m_2,
        A1_h_2,
        A1_rho_2,
        A1_P_2,
        A1_u_2,
        A1_mat_id_2,
    ) = ut.load_particle_data_Fp(
        Fp_init_cond_2, ["pos", "vel", "m", "h", "rho", "P", "u", "mat_id"]
    )

    # Extract only a fof group
    if do_extract_fof_1:
        fof_id = int(im.A1_preset_1[0][4:])

        # Load the group IDs from the fof-output snapshot file
        Fp_fof = (
            remove_substring(im.Fp_init_cond_1, ".hdf5")
            + "_fof_%s_0000.hdf5" % sim.link_len
        )
        A1_fof_id = ut.load_particle_data_Fp(Fp_fof, "fof_id")

        # Select the particles in the group
        A1_sel_fof = np.in1d(A1_fof_id, fof_id)

        print(
            "(1) Extract fof group %d: %d of %d particles"
            % (fof_id, len(A1_sel_fof), len(A1_m_1))
        )
        A2_pos_1 = A2_pos_1[A1_sel_fof]
        A2_vel_1 = A2_vel_1[A1_sel_fof]
        A1_m_1 = A1_m_1[A1_sel_fof]
        A1_h_1 = A1_h_1[A1_sel_fof]
        A1_rho_1 = A1_rho_1[A1_sel_fof]
        A1_P_1 = A1_P_1[A1_sel_fof]
        A1_u_1 = A1_u_1[A1_sel_fof]
        A1_mat_id_1 = A1_mat_id_1[A1_sel_fof]
    if do_extract_fof_2:
        fof_id = int(im.A1_preset_2[0][4:])

        # Load the group IDs from the fof-output snapshot file
        Fp_fof = (
            remove_substring(im.Fp_init_cond_2, ".hdf5")
            + "_fof_%s_0000.hdf5" % sim.link_len
        )
        A1_fof_id = ut.load_particle_data_Fp(Fp_fof, "fof_id")

        # Select the particles in the group
        A1_sel_fof = np.in1d(A1_fof_id, fof_id)

        print(
            "(2) Extract fof group %d: %d of %d particles"
            % (fof_id, len(A1_sel_fof), len(A1_m_2))
        )
        A2_pos_2 = A2_pos_2[A1_sel_fof]
        A2_vel_2 = A2_vel_2[A1_sel_fof]
        A1_m_2 = A1_m_2[A1_sel_fof]
        A1_h_2 = A1_h_2[A1_sel_fof]
        A1_rho_2 = A1_rho_2[A1_sel_fof]
        A1_P_2 = A1_P_2[A1_sel_fof]
        A1_u_2 = A1_u_2[A1_sel_fof]
        A1_mat_id_2 = A1_mat_id_2[A1_sel_fof]

    # Numbers of particles
    num_picle_1 = len(A1_m_1)
    num_picle_2 = len(A1_m_2)
    print(
        "\nNumber of particles: %d + %d = %d = 10^%.3f"
        % (
            num_picle_1,
            num_picle_2,
            num_picle_1 + num_picle_2,
            np.log10(num_picle_1 + num_picle_2),
        )
    )

    # First tweak each set to its own centre of mass and momentum coordinates
    A2_pos_1 -= centre_of_mass(A1_m_1, A2_pos_1)
    A2_vel_1 -= centre_of_mass(A1_m_1, A2_vel_1)
    A2_pos_2 -= centre_of_mass(A1_m_2, A2_pos_2)
    A2_vel_2 -= centre_of_mass(A1_m_2, A2_vel_2)

    # Edit particle data if required
    print("")
    if im.A1_preset_1 is not None:
        for preset in im.A1_preset_1:
            print('Editing (1) with "%s" ' % preset)

            if preset == "v=0":
                A2_vel_1.fill(0)

            elif "rotate=" in preset:
                ax = preset[7]
                angle = float(preset[8:]) * deg_to_rad
                A2_pos_1_old = deepcopy(A2_pos_1)
                A2_vel_1_old = deepcopy(A2_vel_1)

                if ax == "x":
                    A2_pos_1[:, 1] = (
                        np.cos(angle) * A2_pos_1_old[:, 1]
                        - np.sin(angle) * A2_pos_1_old[:, 2]
                    )
                    A2_pos_1[:, 2] = (
                        np.sin(angle) * A2_pos_1_old[:, 1]
                        + np.cos(angle) * A2_pos_1_old[:, 2]
                    )
                    A2_vel_1[:, 1] = (
                        np.cos(angle) * A2_vel_1_old[:, 1]
                        - np.sin(angle) * A2_vel_1_old[:, 2]
                    )
                    A2_vel_1[:, 2] = (
                        np.sin(angle) * A2_vel_1_old[:, 1]
                        + np.cos(angle) * A2_vel_1_old[:, 2]
                    )
                elif ax == "y":
                    A2_pos_1[:, 0] = (
                        np.cos(angle) * A2_pos_1_old[:, 0]
                        + np.sin(angle) * A2_pos_1_old[:, 2]
                    )
                    A2_pos_1[:, 2] = (
                        -np.sin(angle) * A2_pos_1_old[:, 0]
                        + np.cos(angle) * A2_pos_1_old[:, 2]
                    )
                    A2_vel_1[:, 0] = (
                        np.cos(angle) * A2_vel_1_old[:, 0]
                        + np.sin(angle) * A2_vel_1_old[:, 2]
                    )
                    A2_vel_1[:, 2] = (
                        -np.sin(angle) * A2_vel_1_old[:, 0]
                        + np.cos(angle) * A2_vel_1_old[:, 2]
                    )
                elif ax == "z":
                    A2_pos_1[:, 0] = (
                        np.cos(angle) * A2_pos_1_old[:, 0]
                        - np.sin(angle) * A2_pos_1_old[:, 1]
                    )
                    A2_pos_1[:, 1] = (
                        np.sin(angle) * A2_pos_1_old[:, 0]
                        + np.cos(angle) * A2_pos_1_old[:, 1]
                    )
                    A2_vel_1[:, 0] = (
                        np.cos(angle) * A2_vel_1_old[:, 0]
                        - np.sin(angle) * A2_vel_1_old[:, 1]
                    )
                    A2_vel_1[:, 1] = (
                        np.sin(angle) * A2_vel_1_old[:, 0]
                        + np.cos(angle) * A2_vel_1_old[:, 1]
                    )

            elif "spin=" in preset:
                ax = preset[5:]
                A2_pos_1_old = deepcopy(A2_pos_1)
                A2_vel_1_old = deepcopy(A2_vel_1)

                if ax == "x":
                    A2_pos_1[:, 0] = A2_pos_1_old[:, 2]
                    A2_pos_1[:, 2] = -A2_pos_1_old[:, 0]
                    A2_vel_1[:, 0] = A2_vel_1_old[:, 2]
                    A2_vel_1[:, 2] = -A2_vel_1_old[:, 0]
                elif ax == "-x":
                    A2_pos_1[:, 0] = -A2_pos_1_old[:, 2]
                    A2_pos_1[:, 2] = A2_pos_1_old[:, 0]
                    A2_vel_1[:, 0] = -A2_vel_1_old[:, 2]
                    A2_vel_1[:, 2] = A2_vel_1_old[:, 0]
                elif ax == "y":
                    A2_pos_1[:, 1] = A2_pos_1_old[:, 2]
                    A2_pos_1[:, 2] = -A2_pos_1_old[:, 1]
                    A2_vel_1[:, 1] = A2_vel_1_old[:, 2]
                    A2_vel_1[:, 2] = -A2_vel_1_old[:, 1]
                elif ax == "-y":
                    A2_pos_1[:, 1] = -A2_pos_1_old[:, 2]
                    A2_pos_1[:, 2] = A2_pos_1_old[:, 1]
                    A2_vel_1[:, 1] = -A2_vel_1_old[:, 2]
                    A2_vel_1[:, 2] = A2_vel_1_old[:, 1]
                elif ax == "z":
                    pass
                elif ax == "-z":
                    A2_vel_1 *= -1

            elif preset == "xp":
                print("Not using these particles! External potential instead.")

            else:
                raise Exception("Invalid preset", preset)

    if im.A1_preset_2 is not None:
        for preset in im.A1_preset_2:
            print('Editing (2) with "%s" ' % preset)
            if preset == "v=0":
                A2_vel_2.fill(0)

            elif "rotate=" in preset:
                ax = preset[7]
                angle = float(preset[8:]) * deg_to_rad
                A2_pos_2_old = deepcopy(A2_pos_2)
                A2_vel_2_old = deepcopy(A2_vel_2)

                if ax == "x":
                    A2_pos_2[:, 1] = (
                        np.cos(angle) * A2_pos_2_old[:, 1]
                        - np.sin(angle) * A2_pos_2_old[:, 2]
                    )
                    A2_pos_2[:, 2] = (
                        np.sin(angle) * A2_pos_2_old[:, 1]
                        + np.cos(angle) * A2_pos_2_old[:, 2]
                    )
                    A2_vel_2[:, 1] = (
                        np.cos(angle) * A2_vel_2_old[:, 1]
                        - np.sin(angle) * A2_vel_2_old[:, 2]
                    )
                    A2_vel_2[:, 2] = (
                        np.sin(angle) * A2_vel_2_old[:, 1]
                        + np.cos(angle) * A2_vel_2_old[:, 2]
                    )
                elif ax == "y":
                    A2_pos_2[:, 0] = (
                        np.cos(angle) * A2_pos_2_old[:, 0]
                        + np.sin(angle) * A2_pos_2_old[:, 2]
                    )
                    A2_pos_2[:, 2] = (
                        -np.sin(angle) * A2_pos_2_old[:, 0]
                        + np.cos(angle) * A2_pos_2_old[:, 2]
                    )
                    A2_vel_2[:, 0] = (
                        np.cos(angle) * A2_vel_2_old[:, 0]
                        + np.sin(angle) * A2_vel_2_old[:, 2]
                    )
                    A2_vel_2[:, 2] = (
                        -np.sin(angle) * A2_vel_2_old[:, 0]
                        + np.cos(angle) * A2_vel_2_old[:, 2]
                    )
                elif ax == "z":
                    A2_pos_2[:, 0] = (
                        np.cos(angle) * A2_pos_2_old[:, 0]
                        - np.sin(angle) * A2_pos_2_old[:, 1]
                    )
                    A2_pos_2[:, 1] = (
                        np.sin(angle) * A2_pos_2_old[:, 0]
                        + np.cos(angle) * A2_pos_2_old[:, 1]
                    )
                    A2_vel_2[:, 0] = (
                        np.cos(angle) * A2_vel_2_old[:, 0]
                        - np.sin(angle) * A2_vel_2_old[:, 1]
                    )
                    A2_vel_2[:, 1] = (
                        np.sin(angle) * A2_vel_2_old[:, 0]
                        + np.cos(angle) * A2_vel_2_old[:, 1]
                    )

            elif "spin=" in preset:
                ax = preset[5:]
                A2_pos_2_old = deepcopy(A2_pos_2)
                A2_vel_2_old = deepcopy(A2_vel_2)

                if ax == "x":
                    A2_pos_2[:, 0] = A2_pos_2_old[:, 2]
                    A2_pos_2[:, 2] = -A2_pos_2_old[:, 0]
                    A2_vel_2[:, 0] = A2_vel_2_old[:, 2]
                    A2_vel_2[:, 2] = -A2_vel_2_old[:, 0]
                elif ax == "-x":
                    A2_pos_2[:, 0] = -A2_pos_2_old[:, 2]
                    A2_pos_2[:, 2] = A2_pos_2_old[:, 0]
                    A2_vel_2[:, 0] = -A2_vel_2_old[:, 2]
                    A2_vel_2[:, 2] = A2_vel_2_old[:, 0]
                elif ax == "y":
                    A2_pos_2[:, 1] = A2_pos_2_old[:, 2]
                    A2_pos_2[:, 2] = -A2_pos_2_old[:, 1]
                    A2_vel_2[:, 1] = A2_vel_2_old[:, 2]
                    A2_vel_2[:, 2] = -A2_vel_2_old[:, 1]
                elif ax == "-y":
                    A2_pos_2[:, 1] = -A2_pos_2_old[:, 2]
                    A2_pos_2[:, 2] = A2_pos_2_old[:, 1]
                    A2_vel_2[:, 1] = -A2_vel_2_old[:, 2]
                    A2_vel_2[:, 2] = A2_vel_2_old[:, 1]
                elif ax == "z":
                    pass
                elif ax == "-z":
                    A2_vel_2 *= -1

            else:
                raise Exception("Invalid preset", preset)

    # Centre of mass and/or momentum
    if im.is_centre_mass and im.is_centre_mom:
        com = "(centre of mass and momentum)"
    elif im.is_centre_mass:
        com = "(centre of mass)"
    elif im.is_centre_mom:
        com = "(centre of momentum)"
    else:
        com = ""

    # Print scenario parameters
    print("")
    r_unit = sim.po.Di_param_unit_label["r"].unit
    R_print = sim.po.Di_param_unit_label["r"].unit_print
    r_0 = root_sum_sq(im.A1_pos_2_in - im.A1_pos_1_in)
    if im.t_c is not None:
        print("b, B = %.5g, %.5g deg.  " % (im.b, im.B * rad_to_deg), end="")
        print("v_c = %.5g m s^-1 = %.5g v_esc." % (im.v_c, im.v_c_esc))
        print("r_0 = %.5g m, %.5g %s.  " % (r_0, r_0 / r_unit, R_print), end="")
        print("t_c = %.5g s, %.5g h. \n" % (im.t_c, im.t_c * s_to_hour))
    elif im.t_q is not None:
        print("q = %.5g m, %.5g %s.  " % (im.q, im.q / r_unit, R_print), end="")
        if im.v_inf is not None:
            print("v_inf = %.5g m s^-1." % im.v_inf)
        if im.v_q is not None:
            print("v_q = %.5g m s^-1." % im.v_q)
        if im.e is not None:
            print("e = %.5g." % im.e)
        print("r_0 = %.5g m, %.5g %s.  " % (r_0, r_0 / r_unit, R_print), end="")
        print("t_q = %.5g s, %.5g h. \n" % (im.t_q, im.t_q * s_to_hour))

    # Input positions and velocities
    x_1, y_1, z_1 = im.A1_pos_1_in
    v_x_1, v_y_1, v_z_1 = im.A1_vel_1_in
    x_2, y_2, z_2 = im.A1_pos_2_in
    v_x_2, v_y_2, v_z_2 = im.A1_vel_2_in
    print(
        "Initial positions and velocities (input):\n"
        "  (1) [%.5g,  %.5g,  %.5g]  %s \n"
        "      [%.5g,  %.5g,  %.5g]  m s^-1 \n"
        "  (2) [%.5g,  %.5g,  %.5g]  %s \n"
        "      [%.5g,  %.5g,  %.5g]  m s^-1 \n"
        % (
            x_1 / r_unit,
            y_1 / r_unit,
            z_1 / r_unit,
            R_print,
            v_x_1,
            v_y_1,
            v_z_1,
            x_2 / r_unit,
            y_2 / r_unit,
            z_2 / r_unit,
            R_print,
            v_x_2,
            v_y_2,
            v_z_2,
        )
    )

    # Transformed positions and velocities
    if com != "":
        x_1, y_1, z_1 = im.A1_pos_1
        v_x_1, v_y_1, v_z_1 = im.A1_vel_1
        x_2, y_2, z_2 = im.A1_pos_2
        v_x_2, v_y_2, v_z_2 = im.A1_vel_2

        print(
            "Initial positions and velocities %s:\n"
            "  (1) [%.5g,  %.5g,  %.5g]  %s \n"
            "      [%.5g,  %.5g,  %.5g]  m s^-1 \n"
            "  (2) [%.5g,  %.5g,  %.5g]  %s \n"
            "      [%.5g,  %.5g,  %.5g]  m s^-1 \n"
            % (
                com,
                x_1 / r_unit,
                y_1 / r_unit,
                z_1 / r_unit,
                R_print,
                v_x_1,
                v_y_1,
                v_z_1,
                x_2 / r_unit,
                y_2 / r_unit,
                z_2 / r_unit,
                R_print,
                v_x_2,
                v_y_2,
                v_z_2,
            )
        )

    # Shift particle data
    A2_pos_1[:] += im.A1_pos_1
    A2_vel_1[:] += im.A1_vel_1
    A2_pos_2[:] += im.A1_pos_2
    A2_vel_2[:] += im.A1_vel_2

    # Combined particle data
    if "xp" in im.A1_preset_1:
        A2_pos = A2_pos_2
        A2_vel = A2_vel_2
        A1_m = A1_m_2
        A1_h = A1_h_2
        A1_rho = A1_rho_2
        A1_P = A1_P_2
        A1_u = A1_u_2
        A1_mat_id = A1_mat_id_2
    else:
        A2_pos = np.append(A2_pos_1, A2_pos_2, axis=0)
        A2_vel = np.append(A2_vel_1, A2_vel_2, axis=0)
        A1_m = np.append(A1_m_1, A1_m_2)
        A1_h = np.append(A1_h_1, A1_h_2)
        A1_rho = np.append(A1_rho_1, A1_rho_2)
        A1_P = np.append(A1_P_1, A1_P_2)
        A1_u = np.append(A1_u_1, A1_u_2)
        A1_mat_id = np.append(A1_mat_id_1, A1_mat_id_2)
    A1_id = np.arange(len(A2_pos))

    # Save
    Fp_save = sim.Fp_init_cond
    print('Saving to "%s"' % Fp_save[-64:], flush=True)
    with h5py.File(Fp_save, "w") as f:
        woma.save_particle_data(
            f=f,
            A2_pos=A2_pos,
            A2_vel=A2_vel,
            A1_m=A1_m,
            A1_h=A1_h,
            A1_rho=A1_rho,
            A1_P=A1_P,
            A1_u=A1_u,
            A1_mat_id=A1_mat_id,
            A1_id=A1_id,
            boxsize=sim.boxsize,
            file_to_SI=sim.file_to_SI,
        )


def conv_snap_to_init_cond(sim, time):
    """Load a snapshot and save it as initial conditions.

    May be necessary e.g. to run SWIFT's stand-alone FoF when the output
    snapshot files don't include masses and material IDs.

    Parameters
    ----------
    sim : Simulation
        An object with all the information about a simulation.

    time : int
        The snapshot time to load and save. Shortcut -1 for the final snapshot.
    """
    # Set the object if the name was provided instead
    if not isinstance(sim, Simulation):
        sim = Di_sim[sim]
    sim.init()

    # Load the required particle data (including masses etc from the initial conditions)
    (
        A1_id,
        A1_mat_id,
        A1_m,
        A2_pos,
        A2_vel,
        A1_h,
        A1_rho,
        A1_P,
        A1_u,
    ) = sim.load_snapshot_data(
        time, ["id", "mat_id", "m", "pos", "vel", "h", "rho", "P", "u"]
    )

    # Save
    Fp_save = sim.Fp_snap_from_time(time)[:-5] + "_init.hdf5"
    print('Saving to "%s"' % Fp_save[-64:], flush=True)
    with h5py.File(Fp_save, "w") as f:
        woma.save_particle_data(
            f=f,
            A2_pos=A2_pos,
            A2_vel=A2_vel,
            A1_m=A1_m,
            A1_h=A1_h,
            A1_rho=A1_rho,
            A1_P=A1_P,
            A1_u=A1_u,
            A1_mat_id=A1_mat_id,
            A1_id=A1_id,
            boxsize=sim.boxsize,
            file_to_SI=sim.file_to_SI,
        )


def reb_gen_init_cond(rsim):
    """Generate initial conditions and config file for a rebound simulation.

    See utilities.py reb_write_config(), reb_write_init_cond().

    Parameters
    ----------
    rsim : RebSim
        An object with all the information about the REBOUND simulation.
    """
    # Set the object if the name was provided instead
    if not isinstance(rsim, RebSim):
        rsim = Di_reb_sim[rsim]
    rsim.init()

    sim = rsim.sim_in
    sim.init()

    # Set initial positions and velocities of directly provided particles
    for i_picle, p in enumerate(rsim.A1_picle):
        assert p.m is not None
        assert p.R is not None

        # By input position and velocity
        if p.A1_pos is not None:
            assert p.A1_vel is not None
        # By input orbital elements
        elif p.a is not None:
            assert p.e is not None
            assert p.nu is not None

            # Assume particle 0 is the primary
            M_p = rsim.A1_picle[0].m

            # Compute the orbit
            o = Orbit(
                M_p=M_p,
                m=p.m,
                a=p.a,
                e=p.e,
                i=p.i,
                Omega=p.Omega,
                pomega=p.pomega,
                omega=p.omega,
            )

            # Position and velocity at this true anomaly
            p.A1_pos, p.A1_vel = o.pos_vel_from_nu(p.nu)

            # Coordinates relative to the primary
            p.A1_pos += rsim.A2_pos[0]
            p.A1_vel += rsim.A2_vel[0]

        # Set arrays
        rsim.A2_pos[i_picle] = p.A1_pos
        rsim.A2_vel[i_picle] = p.A1_vel

    # Load SPH FoF for additional particles
    if (
        rsim.num_fof is not None
        or rsim.m_fof_min is not None
        or rsim.bnd_orb_Q_max is not None
    ):
        # Load the snapshot particle data
        A1_m, A2_pos, A2_vel, A1_rho, A1_mat_id, A1_fof_id = sim.load_snapshot_data(
            rsim.time_in, ["m", "pos", "vel", "rho", "mat_id", "fof_id"]
        )

        # Reorient by input angle(s)
        for angle, ax in zip(
            [rsim.sim_in_angle, rsim.sim_in_angle_2],
            [rsim.sim_in_rot_ax, rsim.sim_in_rot_ax_2],
        ):
            if angle is not None and angle != 0:
                print("Rotating about %s by %d deg" % (ax, angle))
                angle *= deg_to_rad
                A2_pos_old = deepcopy(A2_pos)
                A2_vel_old = deepcopy(A2_vel)

                if ax == "x":
                    A2_pos[:, 1] = (
                        np.cos(angle) * A2_pos_old[:, 1]
                        - np.sin(angle) * A2_pos_old[:, 2]
                    )
                    A2_pos[:, 2] = (
                        np.sin(angle) * A2_pos_old[:, 1]
                        + np.cos(angle) * A2_pos_old[:, 2]
                    )
                    A2_vel[:, 1] = (
                        np.cos(angle) * A2_vel_old[:, 1]
                        - np.sin(angle) * A2_vel_old[:, 2]
                    )
                    A2_vel[:, 2] = (
                        np.sin(angle) * A2_vel_old[:, 1]
                        + np.cos(angle) * A2_vel_old[:, 2]
                    )
                elif ax == "y":
                    A2_pos[:, 0] = (
                        np.cos(angle) * A2_pos_old[:, 0]
                        + np.sin(angle) * A2_pos_old[:, 2]
                    )
                    A2_pos[:, 2] = (
                        -np.sin(angle) * A2_pos_old[:, 0]
                        + np.cos(angle) * A2_pos_old[:, 2]
                    )
                    A2_vel[:, 0] = (
                        np.cos(angle) * A2_vel_old[:, 0]
                        + np.sin(angle) * A2_vel_old[:, 2]
                    )
                    A2_vel[:, 2] = (
                        -np.sin(angle) * A2_vel_old[:, 0]
                        + np.cos(angle) * A2_vel_old[:, 2]
                    )
                elif ax == "z":
                    A2_pos[:, 0] = (
                        np.cos(angle) * A2_pos_old[:, 0]
                        - np.sin(angle) * A2_pos_old[:, 1]
                    )
                    A2_pos[:, 1] = (
                        np.sin(angle) * A2_pos_old[:, 0]
                        + np.cos(angle) * A2_pos_old[:, 1]
                    )
                    A2_vel[:, 0] = (
                        np.cos(angle) * A2_vel_old[:, 0]
                        - np.sin(angle) * A2_vel_old[:, 1]
                    )
                    A2_vel[:, 1] = (
                        np.sin(angle) * A2_vel_old[:, 0]
                        + np.cos(angle) * A2_vel_old[:, 1]
                    )

        # Reorient by central mass's true anomaly, so sim_in_angle is wrt e.g. the Sun-Mars line
        if rsim.A1_picle[rsim.id_cent].nu is not None:
            angle = rsim.A1_picle[rsim.id_cent].nu
            A2_pos_old = deepcopy(A2_pos)
            A2_vel_old = deepcopy(A2_vel)

            # About z axis
            A2_pos[:, 0] = (
                np.cos(angle) * A2_pos_old[:, 0] - np.sin(angle) * A2_pos_old[:, 1]
            )
            A2_pos[:, 1] = (
                np.sin(angle) * A2_pos_old[:, 0] + np.cos(angle) * A2_pos_old[:, 1]
            )
            A2_vel[:, 0] = (
                np.cos(angle) * A2_vel_old[:, 0] - np.sin(angle) * A2_vel_old[:, 1]
            )
            A2_vel[:, 1] = (
                np.sin(angle) * A2_vel_old[:, 0] + np.cos(angle) * A2_vel_old[:, 1]
            )

        # Initialise to add all fof groups (excluding the no-fof "group")
        fof_id_max = np.amax(A1_fof_id[A1_fof_id < ut.fof_id_none])
        A1_fof_id_add = np.arange(fof_id_max + 1)

        # Group properties
        A1_m_fof = np.empty(len(A1_fof_id_add))
        A2_pos_fof = np.empty((len(A1_fof_id_add), 3))
        A2_vel_fof = np.empty((len(A1_fof_id_add), 3))
        A1_R_fof = np.empty(len(A1_fof_id_add))
        for i_fof, fof_id in enumerate(A1_fof_id_add):
            # Select the particles in the group
            A1_sel = np.in1d(A1_fof_id, fof_id)

            # Mass
            A1_m_fof[i_fof] = sum(A1_m[A1_sel])

            # Centre of mass position and velocity
            A2_pos_fof[i_fof] = centre_of_mass(A1_m[A1_sel], A2_pos[A1_sel])
            A2_vel_fof[i_fof] = centre_of_mass(A1_m[A1_sel], A2_vel[A1_sel])

            # Approximate radius
            A1_R_fof[i_fof] = np.amax(
                root_sum_sq((A2_pos[A1_sel] - A2_pos_fof[i_fof]).T)
            )

        # Remove groups that are below a minimum mass
        if rsim.m_fof_min is not None:
            A1_idx_remove = []
            for i_fof, fof_id in enumerate(A1_fof_id_add):
                if A1_m_fof[i_fof] < rsim.m_fof_min:
                    A1_idx_remove.append(i_fof)

            A1_fof_id_add = np.delete(A1_fof_id_add, A1_idx_remove)
            A1_m_fof = np.delete(A1_m_fof, A1_idx_remove)
            A2_pos_fof = np.delete(A2_pos_fof, A1_idx_remove, axis=0)
            A2_vel_fof = np.delete(A2_vel_fof, A1_idx_remove, axis=0)
            A1_R_fof = np.delete(A1_R_fof, A1_idx_remove)

        # Remove groups that are unbound or above a maximum apoapsis (or custom)
        if rsim.bnd_orb_Q_max is not None:
            A1_idx_remove = []
            for i_fof, fof_id in enumerate(A1_fof_id_add):
                # Estimated orbit
                o = Orbit(
                    A1_pos=A2_pos_fof[i_fof],
                    A1_vel=A2_vel_fof[i_fof],
                    M_p=rsim.A1_m[rsim.id_cent],
                    m=A1_m_fof[i_fof],
                )

                if o.e > 1 or o.Q > rsim.bnd_orb_Q_max:
                    A1_idx_remove.append(i_fof)

            A1_fof_id_add = np.delete(A1_fof_id_add, A1_idx_remove)
            A1_m_fof = np.delete(A1_m_fof, A1_idx_remove)
            A2_pos_fof = np.delete(A2_pos_fof, A1_idx_remove, axis=0)
            A2_vel_fof = np.delete(A2_vel_fof, A1_idx_remove, axis=0)
            A1_R_fof = np.delete(A1_R_fof, A1_idx_remove)

        # Update modified number of fof particles
        if rsim.m_fof_min is not None or rsim.bnd_orb_Q_max is not None:
            if rsim.num_fof in [None, 0] or len(A1_fof_id_add) < rsim.num_fof:
                rsim.num_fof = len(A1_fof_id_add)
                print("Final num_fof = %d" % rsim.num_fof)
    else:
        A1_m_fof = []

    # Test masses
    if rsim.test_masses == "none":
        print("No test particles")
    elif rsim.test_masses == "all":
        print("All orbiting particles as test masses")
    elif rsim.test_masses[:2] == "m<":
        test_mass_cutoff = float(rsim.test_masses[2:])
        print("Test particles for input masses < ", test_mass_cutoff)
    num_test = 0

    # Add additional particles from SPH FoF
    if rsim.num_fof is not None and rsim.num_fof > 0:
        # Shift coordinates to around the central mass
        A2_pos_fof += rsim.A2_pos[rsim.id_cent]
        A2_vel_fof += rsim.A2_vel[rsim.id_cent]

        # Ensure mass order
        A1_sort = np.argsort(A1_m_fof)[::-1]
        A1_m_fof = A1_m_fof[A1_sort]
        A1_R_fof = A1_R_fof[A1_sort]
        A2_pos_fof = A2_pos_fof[A1_sort]
        A2_vel_fof = A2_vel_fof[A1_sort]

        # Add to the particle arrays
        for i_fof in range(rsim.num_fof):
            rsim.A1_m = np.append(rsim.A1_m, A1_m_fof[i_fof])
            rsim.A1_R = np.append(rsim.A1_R, A1_R_fof[i_fof])
            rsim.A2_pos = np.append(rsim.A2_pos, [A2_pos_fof[i_fof]], axis=0)
            rsim.A2_vel = np.append(rsim.A2_vel, [A2_vel_fof[i_fof]], axis=0)

            if rsim.test_masses == "all" or (
                rsim.test_masses[:2] == "m<" and A1_m_fof[i_fof] < test_mass_cutoff
            ):
                num_test += 1

    # Printing
    if len(rsim.A1_picle) > 0:
        print("%d provided particles" % len(rsim.A1_picle))
    if len(A1_m_fof) > 0:
        print(
            "%d FoF particles, masses %.2e to %.2e kg"
            % (len(A1_m_fof), np.amin(A1_m_fof), np.amax(A1_m_fof))
        )
    if num_test > 0:
        print(
            "%d test particles, original masses %.2e to %.2e kg"
            % (num_test, rsim.A1_m[-1], rsim.A1_m[-num_test])
        )

    # Scale up particle radii for collision detection
    if rsim.radius_scale not in [None, 1]:
        print("Scale radii by %g" % rsim.radius_scale)
        rsim.A1_R[rsim.id_cent + 1 :] *= rsim.radius_scale

    # Save the config parameters
    num_picle = len(rsim.A1_m)
    os.makedirs("%s/%s" % (rsim.dir_proj, rsim.name), exist_ok=True)
    ur.reb_write_config(
        Fp_save=rsim.Fp_config,
        num_picle=num_picle,
        num_test=num_test,
        file_init_cond="init_cond.txt",
        file_out_stem="output",
        file_restart="restart.bin",
        t_start=0,
        t_end=rsim.t_end + rsim.t_out_step / 2,
        t_out_step=rsim.t_out_step,
        t_wall_restart=rsim.t_wall_restart,
        t_2=rsim.t_2,
        t_out_step_2=rsim.t_out_step_2,
        t_3=rsim.t_3,
        t_out_step_3=rsim.t_out_step_3,
        id_rel_prim=rsim.id_rel_prim,
        rel_distance_max=rsim.rel_distance_max,
        oblate_id=rsim.oblate_id,
        oblate_J2=rsim.oblate_J2,
        oblate_obliquity=rsim.oblate_obliquity,
        collision_mode=rsim.collision_mode,
        coeff_restitution=rsim.coeff_restitution,
        t_collision_delay=rsim.t_collision_delay,
        id_frame_shift=rsim.id_frame_shift,
    )

    # Save the data
    ur.reb_write_init_cond(
        Fp_save=rsim.Fp_init_cond,
        num_picle=num_picle,
        A1_m=rsim.A1_m,
        A1_R=rsim.A1_R,
        A2_pos=rsim.A2_pos,
        A2_vel=rsim.A2_vel,
    )


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR init_cond.py  ====\n")
