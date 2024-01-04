"""
Specific utilities for GIHR projects: Phobos & Deimos disruptive tidal capture.
"""

from jkeger import *
import gihrpy.utilities as ut


# ========
# Plotting functions
# ========
def func_ax_lim_phodei(po, time, A2_pos):
    """Set axis limits for snapshot plots. See e.g. ut.func_ax_lim_expand."""
    # Custom behaviour for plotting post-disruption fof orbits
    if "fof_orbits" in ut.A1_bonus:
        # Swap the main and inset axes' roles, with zoom-out view for main
        A1_ax_lim = np.array([-241, 60, -293, 8])
        tick_base = 50
        tick_base_minor = None

        # Edit inset axes position
        po.A1_po_inset[0].A1_inset_loc = np.array([0.007, 0.35, 0.643, 0.643])

        # Swap init orbit options
        Di_init_orbit_tmp = deepcopy(po.Di_init_orbit)
        po.Di_init_orbit = po.A1_po_inset[0].Di_init_orbit
        po.A1_po_inset[0].Di_init_orbit = Di_init_orbit_tmp

        # Custom edits for thin panel without any inset
        if ut.A1_arg[0] == "Ma_xp_A2000_n70_r16_v00":
            A1_ax_lim = np.array([-67, 47, -297, 9])
            po.fof_id_max = 30
            po.marker_size = 0.4**2
            po.text_ll = None
            po.text_lr = "scale=50"
            po.A1_po_inset = []

        return A1_ax_lim, tick_base, tick_base_minor

    else:
        return ut.func_ax_lim_expand(po, time, A2_pos)


def func_ax_lim_phodei_inset(po, time, A2_pos):
    """Set axis limits for snapshot insets. See e.g. ut.func_ax_lim_expand."""
    # Custom behaviour for plotting post-disruption fof orbits
    if "fof_orbits" in ut.A1_bonus:
        return ut.func_ax_lim_expand(po, time, A2_pos)

    else:
        return ut.func_ax_lim_expand_step(po, time, A2_pos)


def add_panel_text_disr_snaps(po, time, param_c, ax):
    """Add panel labels for the disruption snapshots figure."""
    time = int(time)

    x = 0.5
    y = 0.96
    ha = "center"
    va = "top"
    c = "k"

    if time == 90000 and param_c == "Q_v_inf":
        text = r"\textbf{e}"
        y = 0.98
    elif time == 90000 and param_c == "fof_id":
        x = 0.06
        y = 0.98
        ha = "left"
        text = r"\textbf{f}"
    else:
        return

    text_outlined(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=po.fontsize,
        color=po.text_colour,
        ha=ha,
        va=va,
        ax=ax,
    )


def plot_Roche_Hill_cross_orbits(po, param_x, param_y, ax):
    """Plot reference orbit lines for the Roche limit and Hill sphere."""
    # Roche limit and Hill sphere
    if param_x == "q" and param_y == "Q":
        ax.axvline(po.R_Roche / po.r_unit, c="k", lw=1.7, ls=ls_dot, zorder=-99)
        ax.axhline(po.R_Hill / po.r_unit, c="k", lw=1.7, ls=ls_dot, zorder=-99)

    # Roche limit and Hill sphere peri/apo
    elif param_x == "a" and param_y == "e":
        A1_e = np.linspace(0, 0.9999, 100)
        A1_a_Roche = po.R_Roche / (1 - A1_e)
        A1_a_Hill = po.R_Hill / (1 + A1_e)

        for A1_a, c in [[A1_a_Roche, "k"], [A1_a_Hill, "0.5"]]:
            ax.plot(
                A1_a / po.r_unit,
                A1_e,
                c=c,
                ls=ls_dot,
                lw=1.7,
                zorder=1,
            )


def add_mass_crossing_text(sim, time, ax):
    """Print mass-crossing text."""
    po = sim.po

    oset = sim.load_snapshot_data(time, "o")

    A1_q = oset.A1_q
    A1_Q = oset.A1_Q
    m_1 = sim.Di_misc["m_c_1"]
    m_2 = sim.Di_misc["m_c_2"]
    m_tot = m_1 + m_2

    # Mass in each region
    m_q_lt_Roche = sum(oset.A1_m[np.where(A1_q < po.R_Roche)[0]])
    m_Q_gt_Hill = sum(oset.A1_m[np.where(A1_Q > po.R_Hill)[0]])
    m_betw_Roche_Hill = sum(
        oset.A1_m[np.where((A1_q > po.R_Roche) & (A1_Q < po.R_Hill))[0]]
    )
    m_betw_1p5Roche_0p5Hill = sum(
        oset.A1_m[np.where((A1_q > 1.5 * po.R_Roche) & (A1_Q < 0.5 * po.R_Hill))[0]]
    )

    m_q_lt_Roche_0 = (m_1 if o_0_1.q < po.R_Roche else 0) + (
        m_2 if o_0_2.q < po.R_Roche else 0
    )
    m_Q_gt_Hill_0 = (m_1 if o_0_1.Q > po.R_Hill else 0) + (
        m_2 if o_0_2.Q > po.R_Hill else 0
    )
    m_betw_Roche_Hill_0 = (
        m_1 if (o_0_1.q > po.R_Roche and o_0_1.Q < po.R_Hill) else 0
    ) + (m_2 if (o_0_2.q > po.R_Roche and o_0_2.Q < po.R_Hill) else 0)
    m_betw_1p5Roche_0p5Hill_0 = (
        m_1 if (o_0_1.q > 1.5 * po.R_Roche and o_0_1.Q < 0.5 * po.R_Hill) else 0
    ) + (m_2 if (o_0_2.q > 1.5 * po.R_Roche and o_0_2.Q < 0.5 * po.R_Hill) else 0)

    # Text
    text = ""
    for i, (label, m, m_0) in enumerate(
        [
            [r"$M_{q < {\rm Roche}}$", m_q_lt_Roche, m_q_lt_Roche_0],
            [r"$M_{Q > {\rm Hill}}$", m_Q_gt_Hill, m_Q_gt_Hill_0],
            [
                r"$M_{q > {\rm Ro}, Q < {\rm Hi}}$",
                m_betw_Roche_Hill,
                m_betw_Roche_Hill_0,
            ],
            [
                r"$M_{q > 1.5 {\rm Ro}, Q < 0.5 {\rm Hi}}$",
                m_betw_1p5Roche_0p5Hill,
                m_betw_1p5Roche_0p5Hill_0,
            ],
        ]
    ):
        if i > 0:
            text += "\n"
        text += label + " = "
        if m == 0:
            text += r"$0$ kg"
        else:
            text += r"$%.2f \!\times\! 10^{%d}$ kg" % (
                m / 10 ** np.floor(np.log10(m)),
                np.floor(np.log10(m)),
            )
        text += (
            r"$\; (%.1f$" % (m_0 / m_tot * 100)
            + r"\%"
            + r"$ \rightarrow %.1f$" % (m / m_tot * 100)
            + r"\%)"
        )

    d = 0.02
    txt = ax.text(
        1 - d,
        d,
        text,
        transform=ax.transAxes,
        fontsize=po.fontsize_text - 8,
        color="k",
        ha="right",
        va="bottom",
    )
    txt.set_path_effects([path_effects.withStroke(linewidth=4, foreground="w")])

    # Scenario info
    text = (
        r"$m_1 = %.2f \!\times\! 10^{%d}$ kg"
        "\n"
        r"$m_2 = %.2f \!\times\! 10^{%d}$ kg"
        "\n"
        r"$r_{\rm c} = %.1f$ %s"
        "\n"
        r"$b = %.1f \; (%.1f^\circ)$"
    ) % (
        m_1 / 10 ** np.floor(np.log10(m_1)),
        np.floor(np.log10(m_1)),
        m_2 / 10 ** np.floor(np.log10(m_2)),
        np.floor(np.log10(m_2)),
        sim.Di_misc["r_c"],
        po.Di_param_unit_label["r"].unit_label,
        sim.Di_misc["b_in"],
        np.arcsin(sim.Di_misc["b_in"]) * rad_to_deg,
    )
    txt = ax.text(
        d,
        1 - d,
        text,
        transform=ax.transAxes,
        fontsize=po.fontsize_text - 8,
        color="k",
        ha="left",
        va="top",
    )
    txt.set_path_effects([path_effects.withStroke(linewidth=4, foreground="w")])


def func_particle_params_extra_phodei(sim, time, param_x, param_y, ax):
    """Misc extras for particle parameters plots."""
    if "init_orbits" in ut.A1_bonus:
        plot_Roche_Hill_cross_orbits(sim.po, param_x, param_y, ax)
        add_mass_crossing_text(sim, time, ax)


def func_param_subsets_extra_phodei(sim_set, param_x, param_y, param_c, param_ls, ax):
    """Misc extras for parameter subset plots."""
    # Dones (1991) predictions
    if (
        sim_set.name == "set_Ma_xp_A2000_s___z_n65_r___v__"
        and param_x == "q"
        and param_c == "v_inf"
    ):
        A1_q = np.linspace(1.01, 3, 100)
        for v_inf in [0.0, 0.2, 0.4, 0.6, 0.8]:
            A1_f_capt_D19 = [
                ut.f_capture_D19(
                    q * Ma.R,
                    v_inf * 1e3,
                    Ma.M,
                    sim_set.A1_sim[0].impact.ic_2.R_s,
                    Ma.o.R_Hill,
                )
                for q in A1_q
            ]

            plt.plot(
                A1_q,
                A1_f_capt_D19,
                color=sim_set.Di_param_Di_p_colour["v_inf"][v_inf],
                ls=ls_dot,
                lw=1.7,
                alpha=0.65,
                zorder=-1,
            )

    # Extra legends
    fs_leg = 22
    if sim_set.name == "set_Ma_xp_A2000_s___z_n65_r___v__":
        # Differentiated asteroids
        A1_leg_extra = [
            ax.plot([], [], c="k", ls=ls_dash_dot_dot, label="Differentiated")[0],
        ]
        leg_extra = ax.legend(
            handles=A1_leg_extra,
            prop={"size": fs_leg},
            loc="upper left",
        )
        ax.add_artist(leg_extra)
    if param_y == "m_capt" and "A2000c30" in sim_set.name:
        # Separate total/core-material
        ax.plot([], [], c="k", label="Total")
        ax.plot([], [], c="k", ls=ls_dash, label="Core")
        ax.legend(prop={"size": fs_leg})

    # Axis tweaks
    if param_x == "q":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.2))
        if sim_set.name == "set_Ma_xp_A2000_n65_r___v0_":
            ax.set_xlim(1.05, 2.25)
        elif sim_set.name in [
            "set_Ma_xp_A2000_s___z_n65_r___v00",
            "set_Ma_xp_A2000_s___z_n65_r___v__",
        ]:
            ax.set_xlim(1.05, 2.85)
    if param_y == "m_capt":
        ax.set_ylabel(r"Mass Fraction", labelpad=-25)


# ========
# REBOUND plotting functions
# ========
def func_reb_orbit_evol_extra_phodei(
    rsim, A1_param_y, time_split, time_split_2, A1_ax_t0, A1_ax_t1, A1_ax_t2
):
    """Misc extras for orbit evolution plots."""
    # Tidy axis ticks etc
    if time_split_2 == 100:
        for ax in A1_ax_t0:
            if time_split <= 6:
                ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
            else:
                ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=2))
                ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=1))
        for ax in A1_ax_t1:
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=20))
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=10))
        for ax in A1_ax_t2:
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1000))
            ax.set_xlim(time_split_2, 5000)

    for i_ax, (ax, param_y) in enumerate(zip(A1_ax_t0, A1_param_y)):
        if param_y == "e":
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.2))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.1))
        elif param_y == "i":
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=30))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=10))
        elif param_y == "q":
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=10))
            # ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=2))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=5))


# ========
# Accumulated results data
# ========
def accum_data_phodei_m_capt(sim):
    """Load pre-recorded values if they exist, or calculate and save them.

    For the Mars initially captured mass etc from SPH tidal disruption simulations.

    Parameters
    ----------
    sim : Simulation
        An object with information about a simulation.

    Returns
    -------
    m_bnd_in_Hill : float
        The mass on bound orbits with apoapses within the Hill sphere.

    m_bnd_in_05Hill : float
        The mass on bound orbits with apoapses within half the Hill sphere.

    m_bnd_out_Hill : float
        The mass on bound orbits with apoapses outside the Hill sphere.

    m_unb : float
        The mass on unbound orbits.

    m_bnd_in_Hill_fof_min : float
        The mass on bound orbits with apoapses within the Hill sphere and in
        fof groups with a mass above m_fof_min.

    num_fof_min : int
        The number of fof groups with a mass above m_fof_min.

    num_fof_min_capt : int
        The number of fof groups with a mass above m_fof_min on initially
        captured orbits.

    m_c_bnd_in_Hill, m_c_bnd_in_05Hill, m_c_bnd_out_Hill, m_c_unb : float
        As above but for core material, for differentiated asteroids.
    """
    print('Loading "%s"' % sim.Fp_accum_data[-64:])
    with h5py.File(sim.Fp_accum_data, "a") as f:
        # Load or compute and save the data
        try:
            m_bnd_in_Hill = f[sim.name].attrs["m_bnd_in_Hill"]
            m_bnd_in_05Hill = f[sim.name].attrs["m_bnd_in_0.5Hill"]
            m_bnd_out_Hill = f[sim.name].attrs["m_bnd_out_Hill"]
            m_unb = f[sim.name].attrs["m_unb"]
            m_bnd_in_Hill_fof_min = f[sim.name].attrs["m_bnd_in_Hill_fof_min"]
            num_fof_min = f[sim.name].attrs["num_fof_min"]
            num_fof_min_capt = f[sim.name].attrs["num_fof_min_capt"]
            modified = f[sim.name].attrs["modified"]

            if "A2000c30" in sim.name:
                m_c_bnd_in_Hill = f[sim.name].attrs["m_c_bnd_in_Hill"]
                m_c_bnd_in_05Hill = f[sim.name].attrs["m_c_bnd_in_0.5Hill"]
                m_c_bnd_out_Hill = f[sim.name].attrs["m_c_bnd_out_Hill"]
                m_c_unb = f[sim.name].attrs["m_c_unb"]

            if modified != "2023_08_24":
                raise KeyError

        except KeyError:
            # Load the final snapshot and fof data
            time = -1
            A1_m, A2_pos, A2_vel, oset, A1_mat_id, A1_fof_id = sim.load_snapshot_data(
                time, ["m", "pos", "vel", "o", "mat_id", "fof_id"]
            )
            oset_fof = load_or_compute_oset_fof(time)

            # Orbits of particles not in fof groups
            A1_sel_not_fof = np.where(A1_fof_id >= ut.fof_id_none)[0]
            if len(A1_sel_not_fof) > 0:
                oset_no_fof = OrbitSet(A1_o=oset.A1_o[A1_sel_not_fof])

            # Calculate the mass subsets
            # Fof groups
            m_bnd_in_Hill = sum(
                oset_fof.A1_m[
                    np.where((oset_fof.A1_e < 1) & (oset_fof.A1_Q < sim.po.R_Hill))[0]
                ]
            )
            m_bnd_in_05Hill = sum(
                oset_fof.A1_m[
                    np.where(
                        (oset_fof.A1_e < 1) & (oset_fof.A1_Q < 0.5 * sim.po.R_Hill)
                    )[0]
                ]
            )
            m_bnd_out_Hill = sum(
                oset_fof.A1_m[np.where(oset_fof.A1_Q >= sim.po.R_Hill)[0]]
            )
            m_unb = sum(oset_fof.A1_m[np.where(oset_fof.A1_e >= 1)[0]])
            m_bnd_in_Hill_fof_min = sum(
                oset_fof.A1_m[
                    np.where(
                        (oset_fof.A1_e < 1)
                        & (oset_fof.A1_Q < sim.po.R_Hill)
                        & (oset_fof.A1_m > 3e15)
                    )[0]
                ]
            )
            num_fof_min = len(np.where(oset_fof.A1_m > 3e15)[0])
            num_fof_min_capt = len(
                np.where(
                    (oset_fof.A1_e < 1)
                    & (oset_fof.A1_Q < sim.po.R_Hill)
                    & (oset_fof.A1_m > 3e15)
                )[0]
            )

            # Non-fof particles
            if len(A1_sel_not_fof) > 0:
                m_bnd_in_Hill += sum(
                    oset_no_fof.A1_m[
                        np.where(
                            (oset_no_fof.A1_e < 1) & (oset_no_fof.A1_Q < sim.po.R_Hill)
                        )[0]
                    ]
                )
                m_bnd_in_05Hill += sum(
                    oset_no_fof.A1_m[
                        np.where(
                            (oset_no_fof.A1_e < 1)
                            & (oset_no_fof.A1_Q < 0.5 * sim.po.R_Hill)
                        )[0]
                    ]
                )
                m_bnd_out_Hill += sum(
                    oset_no_fof.A1_m[np.where(oset_no_fof.A1_Q >= sim.po.R_Hill)[0]]
                )
                m_unb += sum(oset_no_fof.A1_m[np.where(oset_no_fof.A1_e >= 1)[0]])

            # Core material separately
            if "A2000c30" in sim.name:
                A1_f_c_fof = sim.load_fof_data(time, "f_c")
                A1_m_c_fof = oset_fof.A1_m * A1_f_c_fof

                # Fof groups
                m_c_bnd_in_Hill = sum(
                    A1_m_c_fof[
                        np.where((oset_fof.A1_e < 1) & (oset_fof.A1_Q < sim.po.R_Hill))[
                            0
                        ]
                    ]
                )
                m_c_bnd_in_05Hill = sum(
                    A1_m_c_fof[
                        np.where(
                            (oset_fof.A1_e < 1) & (oset_fof.A1_Q < 0.5 * sim.po.R_Hill)
                        )[0]
                    ]
                )
                m_c_bnd_out_Hill = sum(
                    A1_m_c_fof[np.where(oset_fof.A1_Q >= sim.po.R_Hill)[0]]
                )
                m_c_unb = sum(A1_m_c_fof[np.where(oset_fof.A1_e >= 1)[0]])

                # Non-fof particles
                if len(A1_sel_not_fof) > 0:
                    A1_sel_c = A1_mat_id % id_body == sim.mat_id_c
                    m_c_bnd_in_Hill += sum(
                        A1_m[
                            A1_sel_c & A1_sel_not_fof & oset.A1_e
                            < 1 & oset.A1_Q
                            < sim.po.R_Hill
                        ]
                    )
                    m_c_bnd_in_05Hill += sum(
                        A1_m[
                            A1_sel_c & A1_sel_not_fof & oset.A1_e
                            < 1 & oset.A1_Q
                            < 0.5 * sim.po.R_Hill
                        ]
                    )
                    m_c_bnd_out_Hill += sum(
                        A1_m[A1_sel_c & A1_sel_not_fof & oset.A1_Q >= sim.po.R_Hill]
                    )
                    m_c_unb += sum(A1_m[A1_sel_c & A1_sel_not_fof & oset.A1_e >= 1])

            # Create/select the group and write the data
            grp = f.require_group("/%s" % sim.name)
            print('Writing "%s/%s" etc...' % (sim.name, "m_bnd_in_Hill"))
            grp.attrs["m_bnd_in_Hill"] = m_bnd_in_Hill
            grp.attrs["m_bnd_in_0.5Hill"] = m_bnd_in_05Hill
            grp.attrs["m_bnd_out_Hill"] = m_bnd_out_Hill
            grp.attrs["m_unb"] = m_unb
            grp.attrs["m_bnd_in_Hill_fof_min"] = m_bnd_in_Hill_fof_min
            grp.attrs["num_fof_min"] = num_fof_min
            grp.attrs["num_fof_min_capt"] = num_fof_min_capt
            grp.attrs["modified"] = "2023_08_24"
            if "A2000c30" in sim.name:
                grp.attrs["m_c_bnd_in_Hill"] = m_c_bnd_in_Hill
                grp.attrs["m_c_bnd_in_0.5Hill"] = m_c_bnd_in_05Hill
                grp.attrs["m_c_bnd_out_Hill"] = m_c_bnd_out_Hill
                grp.attrs["m_c_unb"] = m_c_unb

        # Return
        if "A2000c30" in sim.name:
            return (
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
            )
        else:
            return (
                m_bnd_in_Hill,
                m_bnd_in_05Hill,
                m_bnd_out_Hill,
                m_unb,
                m_bnd_in_Hill_fof_min,
                num_fof_min,
                num_fof_min_capt,
            )


def accum_data_phodei_reb_m(rsim):
    """Load pre-recorded values if they exist, or calculate and save them.

    For the Mars final surviving mass etc from REBOUND integration simulations.

    Parameters
    ----------
    rsim : RebSim
        The rebound simulation object.

    Returns
    -------
    m_tot_in : float
        The total input mass.

    m_rm_close : float
        The mass removed by colliding with Mars.

    m_rm_far : float
        The mass removed by escaping far beyond Mars's Hill sphere.

    m_bump : float
        The mass that experienced a low-energy, non-disruptive collision.

    m_disr : float
        The mass that suffered a disruptive collision.

    m_no_col : int
        The mass that had no collisions.

    m_surv : int
        The mass still in orbit.

    m_col : float
        The mass that survived in orbit or suffered a disruptive collision.
    """
    print('Loading "%s"' % rsim.Fp_accum_data[-64:])
    with h5py.File(rsim.Fp_accum_data, "a") as f:
        # Load or compute and save the data
        try:
            m_tot_in = f[rsim.name].attrs["m_tot_in"]
            m_rm_close = f[rsim.name].attrs["m_rm_close"]
            m_rm_far = f[rsim.name].attrs["m_rm_far"]
            m_bump = f[rsim.name].attrs["m_bump"]
            m_disr = f[rsim.name].attrs["m_disr"]
            m_no_col = f[rsim.name].attrs["m_no_col"]
            m_surv = f[rsim.name].attrs["m_surv"]
            m_col = f[rsim.name].attrs["m_col"]

            if f[rsim.name].attrs["modified"] != "2023_09_29":
                raise KeyError

        except KeyError:
            # Load output data
            (
                A1_time,
                A1_m,
                A1_R,
                A3_pos,
                A3_vel,
                A1_idx_last,
                A1_fate,
            ) = rsim.reb_load_data()
            num_orb = len(A1_m) - rsim.id_orb_start
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

            # Initial mass of (assume bound) orbiting particles
            A1_m_orb = A1_m[rsim.id_orb_start :]
            m_tot_in = sum(A1_m_orb)
            assert rsim.bnd_orb_Q_max is not None

            # Surviving particles
            A1_sel_surv = np.where(A1_idx_last[rsim.id_orb_start :] == -1)[0]

            # Particles that have a disruptive collision
            if A1_id_1_c is None or len(A1_id_1_c) == 0:
                A1_sel_col = []
                A1_sel_no_col = np.arange(len(A1_m_orb))
                A1_sel_bump = []
            else:
                # Any collision
                # A1_sel_col = np.unique(np.append(A1_id_1_c, A1_id_2_c)) - rsim.id_orb_start

                # Require a disruptive collision
                A1_m_1_c = np.array([A1_m[id] for id in A1_id_1_c])
                A1_m_2_c = np.array([A1_m[id] for id in A1_id_2_c])
                A1_R_1_c = np.array([A1_R[id] for id in A1_id_1_c])
                A1_R_2_c = np.array([A1_R[id] for id in A1_id_2_c])
                A1_m_tot = A1_m_1_c + A1_m_2_c
                A1_mu = A1_m_1_c * A1_m_2_c / A1_m_tot

                # Collision energy
                A1_Q_c = 0.5 * A1_mu * A1_v_rel_c**2 / A1_m_tot
                A1_Q_star_RD = np.array(
                    [
                        Q_V_star_RD_SL12(
                            M_t=A1_m_2_c[i_c],
                            M_i=A1_m_1_c[i_c],
                            R_t=A1_R_2_c[i_c],
                            R_i=A1_R_1_c[i_c],
                            beta=45 * deg_to_rad,
                            c_star=5,
                            mu_=0.36,
                        )[0]
                        for i_c in range(len(A1_time_c))
                    ]
                )

                # Highest collision energy of each particle
                A1_Q_c_max = np.zeros_like(A1_m_orb)
                for i_orb in range(len(A1_m_orb)):
                    i_picle = rsim.id_orb_start + i_orb
                    A1_sel_i_orb = np.where(
                        (A1_id_1_c == i_picle) | (A1_id_2_c == i_picle)
                    )[0]

                    if len(A1_sel_i_orb) > 0:
                        A1_Q_c_max[i_orb] = np.amax(
                            (A1_Q_c / A1_Q_star_RD)[A1_sel_i_orb]
                        )

                # Particles with at least one disruptive collision
                A1_sel_col = np.where(A1_Q_c_max > 1)[0]

                # Particles with no collisions
                A1_sel_no_col = np.where(A1_Q_c_max == 0)[0]

                # Particles with only non-disruptive collisions
                A1_sel_bump = np.where((A1_Q_c_max > 0) & (A1_Q_c_max < 1))[0]

            # Surviving and/or colliding particles
            A1_sel_col_surv = np.unique(np.append(A1_sel_col, A1_sel_surv)).astype(int)

            # Masses by fate
            A1_rm_close = np.where(
                A1_fate[rsim.id_orb_start :] == cl.RebSim.reb_fate_rm_close
            )[0]
            A1_rm_far = np.where(
                A1_fate[rsim.id_orb_start :] == cl.RebSim.reb_fate_rm_far
            )[0]

            m_rm_close = sum(A1_m_orb[A1_rm_close])
            m_rm_far = sum(A1_m_orb[A1_rm_far])
            m_bump = sum(A1_m_orb[A1_sel_bump])
            m_disr = sum(A1_m_orb[A1_sel_col])
            m_no_col = sum(A1_m_orb[A1_sel_no_col])
            m_surv = sum(A1_m_orb[A1_sel_surv])
            m_col = sum(A1_m_orb[A1_sel_col_surv])

            # Create/select the group and write the data
            grp = f.require_group("/%s" % rsim.name)
            print('Writing "%s/%s" etc...' % (rsim.name, "m_surv"))
            grp.attrs["m_tot_in"] = m_tot_in
            grp.attrs["m_rm_close"] = m_rm_close
            grp.attrs["m_rm_far"] = m_rm_far
            grp.attrs["m_bump"] = m_bump
            grp.attrs["m_disr"] = m_disr
            grp.attrs["m_no_col"] = m_no_col
            grp.attrs["m_surv"] = m_surv
            grp.attrs["m_col"] = m_col
            grp.attrs["modified"] = "2023_09_29"

        return m_tot_in, m_rm_close, m_rm_far, m_bump, m_disr, m_no_col, m_surv, m_col


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR utilities_phodei.py  ====\n")
