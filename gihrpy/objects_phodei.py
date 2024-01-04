"""
Objects for GIHR: Phobos & Deimos disruptive tidal capture project.

Note that some blocks are indented in `if True:` statements purely for
convenient block-folding in the code editor.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.utilities_reb as ur
import gihrpy.utilities_phodei as ut_ph
from gihrpy.classes import PlotOptions, InitProf, InitCond
from gihrpy.classes_sim import ImpactInitCond, Simulation, SimSet
from gihrpy.classes_reb import RebPicle, RebSim


# ========
# Naming elements
# ========
# Ma        Mars
# AXXYY     Asteroid mass, 10^XX.YY kg  (cZZ core fraction ZZ%)
# sXXYZ     Spin period, XX.Y hours, ang-mom in Z direction (m=minus)
# nXY       Number of particles in the body, 10^X.Y
# xp        External potential (instead of nXY)
# rXY       Periapsis distance, q = X.Y Mars radii (Roche ~2.6)
# vXY       Speed at infinity, v_inf = X.Y km/s

# phiXXX    Angle of encounter periapsis relative to Mars's velocity, XXX deg
# iXX       Inclination of encounter, XX deg
# oXX       Obliquity of Mars, XX deg
# eXYY      Eccentricity of Mars, X.YY
# nuXXX     Initial true anomaly of Mars, XXX deg

# Reference masses
#           (kg)        (Mars)      (Earth)     (Phobos)  (Deimos)  (P + D)
# Mars      6.41e23     1           0.1074      5.99e7    4.33e8    5.27e7
# Phobos    1.07e16     1.66e-8     1.78e-10    1         7.23      0.878
# Deimos    1.48e15     2.30e-9     2.47e-10    0.138     1         0.122
# P + D     1.22e16     1.90e-8     2.04e-9     1.138     8.23      1

# Plotting Options
if True:
    plot_opt_phodei         = PlotOptions(
        dir_save            = ut.dir_plots + "phodei/",
        Di_mat_colour       = {
            "ANEOS_Fe85Si15"    : "#808080",
            "ANEOS_forsterite"  : "#cc3300",
            "ANEOS_Fe85Si15_2"  : "#775533",
            # "ANEOS_forsterite_2": "#eecc00",
            "ANEOS_forsterite_2": "#4400bb",
        },
        Di_misc_colour      = {
            "planet"    : "#aa3300",
            "fof_small" : "#330044",
            "fof_none"  : "0.4",
            "Roche"     : "k",
            "Hill"      : "0.88",
            "o_1"       : "#cc3300",
            "o_2"       : "#4400bb",
            "com"       : "#008800",
        },
        Di_mat_highlt       = {
            "ANEOS_Fe85Si15"    : "#4433aa",
            "ANEOS_forsterite"  : "#7711dd",
            "ANEOS_Fe85Si15_2"  : "#3344aa",
            "ANEOS_forsterite_2": "#1177dd",
        },
        unit_label_r_def    = ParamUnitLabel("r", "Distance", Ma.R, r"$R_{\mars{}\!}$", "R_M", " * Ma.R"),
        unit_label_m_def    = ParamUnitLabel("m", "Mass", 1, "kg", "kg", ""),
        Di_param_unit_label = {
            "v_inf": ParamUnitLabel("v_inf", r"$v_\infty$", 1e3, r"km s$^{-1}$"),
        },
        Di_param_lim        = {
            "r"         : [0.5, 10.5],
            # "a"         : [1, 71],
            "a"         : [2, 300],
            "e"         : [0, 1.005],
            "q"         : [1, 30],
            "Q"         : [0, 400],
            "a_eq"      : [0, 11],
            "m_capt"    : [0.05, 1],
        },
        Di_param_log        = {
            "m_capt"    : True,
        },
        func_cmap_param     = ut.func_cmap_param_def,
        A1_ax_lim           = np.array([-0.92, 1.08, -1.08, 0.92]) * 0.35,
        tick_base           = 0.2,
        func_ax_lim         = ut_ph.func_ax_lim_phodei,
        Di_ax_lim           = {
            "f_extent"      : 5e-3,
            "pad"           : 0.27,
        },
        thick               = -0.0001,
        # thick             = 0,
        # thick             = "None",
        func_marker_size    = ut.func_marker_size_ax_scale,
        Di_marker_size      = {
            "marker_size_ref"   : 1.2**2,
            "ax_size_ref"       : 0.1,
            "pow_size_scale"    : 0.7,
        },
        ref_m_picle         = 10**20 / 10**7,
        alpha               = 0.3,
        Di_init_orbit       = {
            "colour"        : "k",
            "lw"            : 3.5,
            "alpha"         : 0.4,
            "ls_step_on"    : 0.01,
            "ls_step_off"   : 0.21,
        },
        fof_id_max          = 49,
        fof_n_min           = 20,
        R_Roche             = 2.63 * Ma.R,
        R_Hill              = Ma.o.R_Hill,
        text_ul             = "tstamp",
        text_ll             = "scale=auto",
        func_particle_params_extra  = ut_ph.func_particle_params_extra_phodei,
        func_param_subsets_extra    = ut_ph.func_param_subsets_extra_phodei,
    )
    plot_opt_phodei_init    = PlotOptions(
        copy                = plot_opt_phodei,
        unit_label_r_def    = ParamUnitLabel("r", "Distance", 1e3, "kg", "kg", " * 1e3"),
        A1_ax_lim           = np.array([-1, 1, -1, 1]) * 500,
        tick_base           = "None",
        thick               = 0,
        marker_size         = 2**2,
        alpha               = 0.6,
    )
    plot_opt_phodei_reb     = PlotOptions(
        copy            = plot_opt_phodei,
        Di_param_lim        = {
            "e"         : [0, 1],
            "i"         : [0, 83],
            "q"         : [0, 35],
        },
        func_reb_orbit_evol_extra   = ut_ph.func_reb_orbit_evol_extra_phodei,
    )

    # Tweaked options after copies
    plot_opt_phodei.A1_po_inset = [
        PlotOptions(
            copy                = plot_opt_phodei,
            A1_inset_loc        = [0.65, 0.007, 0.343, 0.343],
            A1_ax_lim           = np.array([-18, 7, -20, 5]),
            func_ax_lim         = ut_ph.func_ax_lim_phodei_inset,
            Di_ax_lim           = {
                "f_extent"      : plot_opt_phodei.Di_ax_lim["f_extent"],
                "pad"           : plot_opt_phodei.Di_ax_lim["pad"] + 0.05,
                "step"          : 10,
            },
            Di_marker_size      = {
                "marker_size_ref"   : plot_opt_phodei.Di_marker_size["marker_size_ref"],
                "ax_size_ref"       : plot_opt_phodei.Di_marker_size["ax_size_ref"],
                "pow_size_scale"    : 0.1,
            },
            Di_init_orbit       = {
                "colour"        : "k",
                "lw"            : 1,
                "alpha"         : 0.4,
                "ls"            : "-",
                # "ls_step_on"    : 0.01,
                # "ls_step_off"   : 0.035,
            },
        )
    ]

# Profiles
if True:
    prof_Ma         = InitProf(
        name            = "prof_Ma",
        dir_proj        = ut.dir_gihr + "swift/phodei/",
        po              = plot_opt_phodei_init,
        num_prof        = 10000,
        A1_mat_layer    = ["ANEOS_Fe85Si15", "ANEOS_forsterite"],
        A1_T_rho_type   = ["adiabatic", "adiabatic"],
        R_min           = 0.95 * Ma.R,
        R_max           = 1.05 * Ma.R,
        # A1_M_layer    = np.array([0.24, 0.76]) * Ma.M,
        # P_s           = 1e5,
        # T_s           = 1000,
        A1_M_layer      = np.array([0.24, 0.76]) * Ma.M,
        A1_R_layer      = np.array([0.48174, 0.99791]) * Ma.R,
        A1_idx_layer    = [4827, 9999],
        M               = 1 * Ma.M,
        P_s             = 1e+05,
        T_s             = 1000,
        rho_s           = 3164.4,
        P_1             = 2.1468e+10,
        T_1             = 1106.5,
        rho_1           = 8118.8,
        P_0             = 4.9079e+10,
        T_0             = 1298.5,
        rho_0           = 8897.5,
        I_MR2           = 0.35079,
        Di_po_edit      = {
            "unit_label_r_def"  : ParamUnitLabel("r", "Distance", Ma.R, "$R_{\mars{}\!}$", "R_M", " * Ma.R"),
            "unit_label_m_def"  : ParamUnitLabel("m", "Mass", Ma.M, r"$M_{\mars{}\!}$", "M_M", " * Ma.M"),
        }
    )

    # Masses
    prof_A1800      = InitProf(
        name            = "prof_A1800",
        dir_proj        = ut.dir_gihr + "swift/phodei/",
        po              = plot_opt_phodei_init,
        num_prof        = 10000,
        A1_mat_layer    = ["ANEOS_forsterite"],
        A1_T_rho_type   = ["power=0"],
        R_max           = 500e3,
        # M               = 10**18.00,
        # P_s             = 1e4,
        # T_s             = 500,
        A1_M_layer      = np.array([9.9975e+17]),
        A1_R_layer      = np.array([42.049]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 9.9975e+17,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 2.5623e+06,
        T_0             = 500,
        rho_0           = 3209.7,
        I_MR2           = 0.40004,
    )
    prof_A1850      = InitProf(
        name            = "prof_A1850",
        copy            = "prof_A1800",
        # M               = 10**18.50,
        A1_M_layer      = np.array([3.1618e+18]),
        A1_R_layer      = np.array([61.722]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 3.1618e+18,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 5.5066e+06,
        T_0             = 500,
        rho_0           = 3209.8,
        I_MR2           = 0.40004,
    )
    prof_A1900      = InitProf(
        name            = "prof_A1900",
        copy            = "prof_A1800",
        # M               = 10**19.00,
        A1_M_layer      = np.array([1e+19]),
        A1_R_layer      = np.array([90.599]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 1e+19,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 1.1848e+07,
        T_0             = 500,
        rho_0           = 3209.9,
        I_MR2           = 0.40004,
    )
    prof_A1950      = InitProf(
        name            = "prof_A1950",
        copy            = "prof_A1800",
        # M               = 10**19.50,
        A1_M_layer      = np.array([3.1623e+19]),
        A1_R_layer      = np.array([132.98]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 3.1623e+19,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 2.5518e+07,
        T_0             = 500,
        rho_0           = 3209.9,
        I_MR2           = 0.40004,
    )
    prof_A2000      = InitProf(
        name            = "prof_A2000",
        copy            = "prof_A1800",
        # M               = 10**20.00,
        A1_M_layer      = np.array([9.9981e+19]),
        A1_R_layer      = np.array([195.18]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 9.9981e+19,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 5.4927e+07,
        T_0             = 500,
        rho_0           = 3210,
        I_MR2           = 0.40004,
    )
    prof_A2050      = InitProf(
        name            = "prof_A2050",
        copy            = "prof_A1800",
        # M               = 10**20.50,
        A1_M_layer      = np.array([3.1621e+20]),
        A1_R_layer      = np.array([286.48]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 3.1621e+20,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 1.1848e+08,
        T_0             = 500,
        rho_0           = 3211.9,
        I_MR2           = 0.40002,
    )
    prof_A2100      = InitProf(
        name            = "prof_A2100",
        copy            = "prof_A1800",
        # M               = 10**21.00,
        A1_M_layer      = np.array([9.9995e+20]),
        A1_R_layer      = np.array([420.44]) * 1e3,
        A1_idx_layer    = [9999],
        M               = 9.9995e+20,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_0             = 2.5529e+08,
        T_0             = 500,
        rho_0           = 3215.4,
        I_MR2           = 0.39995,
    )

    # Spins  (~[1, 3/4, 1/2, 1/4, 1/8] L_max)  (encounter orbit L ~ 2e29)
    prof_A2000_s030 = InitProf(
        name            = "prof_A2000_s030",
        planet          = prof_A2000,
        period          = 3.0,
        A1_M_layer      = np.array([9.9509e+19]),
        A1_R_layer      = np.array([230.9]) * 1e3,
        A1_Z_layer      = np.array([138.81]) * 1e3,
        A1_idx_layer    = [1774],
        M               = 9.9509e+19,
        P_s             = 44493,
        T_s             = 500,
        rho_s           = 3209.4,
        P_0             = 3.6209e+07,
        T_0             = 500,
        rho_0           = 3209.9,
        I_MR2           = 0.4,
        L               = 1.2346e+27, # 3.53e-08 L_EM
        L_max           = 1.23e27,
    )
    prof_A2000_s036 = InitProf(
        name            = "prof_A2000_s036",
        planet          = prof_A2000,
        period          = 3.6,
        A1_M_layer      = np.array([1.0016e+20]),
        A1_R_layer      = np.array([216.72]) * 1e3,
        A1_Z_layer      = np.array([158.61]) * 1e3,
        A1_idx_layer    = [1665],
        M               = 1.0016e+20,
        P_s             = 21556,
        T_s             = 500,
        rho_s           = 3209.4,
        P_0             = 4.3342e+07,
        T_0             = 500,
        rho_0           = 3210,
        I_MR2           = 0.4,
        L               = 9.1225e+26, # 2.61e-08 L_EM
        L_max           = 1.23e27,
    )
    prof_A2000_s047 = InitProf(
        name            = "prof_A2000_s047",
        planet          = prof_A2000,
        period          = 4.7,
        A1_M_layer      = np.array([1.0007e+20]),
        A1_R_layer      = np.array([205]) * 1e3,
        A1_Z_layer      = np.array([177.1]) * 1e3,
        A1_idx_layer    = [1575],
        M               = 1.0007e+20,
        P_s             = 55196,
        T_s             = 500,
        rho_s           = 3209.5,
        P_0             = 4.8847e+07,
        T_0             = 500,
        rho_0           = 3210,
        I_MR2           = 0.4,
        L               = 6.2468e+26, # 1.78e-08 L_EM
        L_max           = 1.23e27,
    )
    prof_A2000_s086 = InitProf(
        name            = "prof_A2000_s086",
        planet          = prof_A2000,
        period          = 8.6,
        A1_M_layer      = np.array([1e+20]),
        A1_R_layer      = np.array([196.8]) * 1e3,
        A1_Z_layer      = np.array([192.04]) * 1e3,
        A1_idx_layer    = [1512],
        M               = 1e+20,
        P_s             = 31627,
        T_s             = 500,
        rho_s           = 3209.4,
        P_0             = 5.3147e+07,
        T_0             = 500,
        rho_0           = 3210,
        I_MR2           = 0.4,
        L               = 3.1442e+26, # 8.98e-09 L_EM
        L_max           = 1.23e27,
    )
    prof_A2000_s170 = InitProf(
        name            = "prof_A2000_s170",
        planet          = prof_A2000,
        period          = 17.0,
        A1_M_layer      = np.array([9.9851e+19]),
        A1_R_layer      = np.array([195.5]) * 1e3,
        A1_Z_layer      = np.array([194.3]) * 1e3,
        A1_idx_layer    = [1502],
        M               = 9.9851e+19,
        P_s             = 33889,
        T_s             = 500,
        rho_s           = 3209.4,
        P_0             = 5.4423e+07,
        T_0             = 500,
        rho_0           = 3210,
        I_MR2           = 0.4,
        L               = 1.5672e+26, # 4.48e-09 L_EM
        L_max           = 1.23e27,
    )

    # Differentiated
    prof_A2000c30   = InitProf(
        name            = "prof_A2000c30",
        copy            = "prof_A2000",
        A1_mat_layer    = ["ANEOS_Fe85Si15", "ANEOS_forsterite"],
        A1_T_rho_type   = ["power=0", "power=0"],
        # A1_M_layer      = np.array([0.3, 0.7]) * 10**20.00,
        R_min           = 170e3,
        R_max           = 190e3,
        A1_M_layer      = np.array([2.999e+19, 7.0002e+19]),
        A1_R_layer      = np.array([98.726, 183.39]) * 1e3,
        A1_idx_layer    = [5383, 9999],
        M               = 9.9992e+19,
        P_s             = 10000,
        T_s             = 500,
        rho_s           = 3209.3,
        P_1             = 5.1517e+07,
        T_1             = 500,
        rho_1           = 7438.1,
        P_0             = 1.2716e+08,
        T_0             = 500,
        rho_0           = 7438.4,
        I_MR2           = 0.35156,
    )

# Initial particles
if True:
    init_Ma_n50     = InitCond(
        name            = "init_Ma_n50",
        ip              = prof_Ma,
        file_to_SI      = Conversions(1e24, 1e6, 1),
        boxsize         = 2 * 5.0 * 1e6,
        num_picle_des   = 10**5,
        num_picle       = 111835,
        M               = 0.999698 * Ma.M,
        R_s             = 0.997815 * Ma.R,
        min_sep         = 86.7 * 1e3,
    )

    # Masses
    init_A1800_n60  = InitCond(
        name            = "init_A1800_n60",
        ip              = prof_A1800,
        file_to_SI      = Conversions(1e24, 1e6, 1),
        boxsize         = 2 * 500e3,
        num_picle_des   = 10**6,
        num_picle       = 1002827,
        M               = 9.99557e+17,
        R_s             = 42.0494 * 1e3,
        min_sep         = 0.679 * 1e3,
    )
    init_A1850_n60  = InitCond(
        name            = "init_A1850_n60",
        copy            = "init_A1800_n60",
        ip              = prof_A1850,
        num_picle_des   = 10**6,
        num_picle       = 1002827,
        M               = 3.1612e+18,
        R_s             = 61.7218 * 1e3,
        min_sep         = 0.997 * 1e3,
    )
    init_A1900_n65  = InitCond(
        name            = "init_A1900_n65",
        copy            = "init_A1800_n60",
        ip              = prof_A1900,
        num_picle_des   = 10**6.5,
        num_picle       = 3227694,
        M               = 9.99816e+18,
        R_s             = 90.5991 * 1e3,
        min_sep         = 0.992 * 1e3,
    )
    init_A1950_n65  = InitCond(
        name            = "init_A1950_n65",
        copy            = "init_A1800_n60",
        ip              = prof_A1950,
        num_picle_des   = 10**6.5,
        num_picle       = 3227694,
        M               = 3.16169e+19,
        R_s             = 132.98 * 1e3,
        min_sep         = 1.46 * 1e3,
    )
    init_A2000_n65  = InitCond(
        name            = "init_A2000_n65",
        copy            = "init_A1800_n60",
        ip              = prof_A2000,
        num_picle_des   = 10**6.5,
        num_picle       = 3227695,
        M               = 9.99627e+19,
        R_s             = 195.175 * 1e3,
        min_sep         = 2.14 * 1e3,
    )
    init_A2050_n70  = InitCond(
        name            = "init_A2050_n70",
        copy            = "init_A1800_n60",
        ip              = prof_A2050,
        num_picle_des   = 10**7.0,
        num_picle       = 10094432,
        M               = 3.16156e+20,
        R_s             = 286.484 * 1e3,
        min_sep         = 2.14 * 1e3,
    )
    init_A2100_n70  = InitCond(
        name            = "init_A2100_n70",
        copy            = "init_A1800_n60",
        ip              = prof_A2100,
        num_picle_des   = 10**7.0,
        num_picle       = 10087874,
        M               = 9.99782e+20,
        R_s             = 420.441 * 1e3,
        min_sep         = 3.15 * 1e3,
    )

    # Spins
    init_A2000_s030_n65 = InitCond(
        name            = "init_A2000_s030_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000_s030,
        num_picle       = 3162277,
        M               = 9.95094e+19,
        R_s             = 230.904 * 1e3,
        Z_s             = 138.812 * 1e3,
        R_s_sph         = 158.418 * 1e3,
        min_sep         = 1.67 * 1e3,
    )
    init_A2000_s036_n65 = InitCond(
        name            = "init_A2000_s036_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000_s036,
        num_picle       = 3162274,
        M               = 1.00161e+20,
        R_s             = 216.717 * 1e3,
        Z_s             = 158.614 * 1e3,
        R_s_sph         = 173.295 * 1e3,
        min_sep         = 1.87 * 1e3,
    )
    init_A2000_s047_n65 = InitCond(
        name            = "init_A2000_s047_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000_s047,
        num_picle       = 3162278,
        M               = 1.0007e+20,
        R_s             = 205.002 * 1e3,
        Z_s             = 177.097 * 1e3,
        R_s_sph         = 183.963 * 1e3,
        min_sep         = 2.14 * 1e3,
    )
    init_A2000_s086_n65 = InitCond(
        name            = "init_A2000_s086_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000_s086,
        num_picle       = 3162278,
        M               = 1.00005e+20,
        R_s             = 196.802 * 1e3,
        Z_s             = 192.037 * 1e3,
        R_s_sph         = 191.91 * 1e3,
        min_sep         = 2.13 * 1e3,
    )
    init_A2000_s170_n65 = InitCond(
        name            = "init_A2000_s170_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000_s170,
        num_picle       = 3162275,
        M               = 9.98513e+19,
        R_s             = 195.501 * 1e3,
        Z_s             = 194.303 * 1e3,
        R_s_sph         = 194.185 * 1e3,
        min_sep         = 2.15 * 1e3,
    )

    # Differentiated
    init_A2000c30_n65   = InitCond(
        name            = "init_A2000c30_n65",
        copy            = "init_A2000_n65",
        ip              = prof_A2000c30,
        num_picle       = 3195625,
        M               = 9.9963e+19,
        R_s             = 183.385 * 1e3,
        min_sep         = 1.62 * 1e3,
    )

    # Resolutions
    init_A2000_n50  = InitCond(
        name            = "init_A2000_n50",
        copy            = "init_A2000_n65",
        num_picle_des   = 10**5.0,
        num_picle       = 108605,
        M               = 9.99557e+19,
        R_s             = 195.175 * 1e3,
        min_sep         = 6.62 * 1e3,
    )
    init_A2000_n55  = InitCond(
        name            = "init_A2000_n55",
        copy            = "init_A2000_n65",
        num_picle_des   = 10**5.5,
        num_picle       = 319379,
        M               = 9.99597e+19,
        R_s             = 195.175 * 1e3,
        min_sep         = 4.62 * 1e3,
    )
    init_A2000_n60  = InitCond(
        name            = "init_A2000_n60",
        copy            = "init_A2000_n65",
        num_picle_des   = 10**6.0,
        num_picle       = 1002827,
        M               = 9.99615e+19,
        R_s             = 195.175 * 1e3,
        min_sep         = 3.15 * 1e3,
    )
    init_A2000_n70  = InitCond(
        name            = "init_A2000_n70",
        copy            = "init_A2000_n65",
        num_picle_des   = 10**7.0,
        num_picle       = 10099881,
        M               = 9.99636e+19,
        R_s             = 195.175 * 1e3,
        min_sep         = 1.46 * 1e3,
    )

# Settling simulations
if True:
    Ma_n50          = Simulation(
        name            = "Ma_n50",
        ic              = init_Ma_n50,
        po              = plot_opt_phodei_init,
        dir_proj        = ut.dir_gihr + "swift/phodei/",
        category        = "phodei_init",
        file_to_SI      = Conversions(1e24, 1e6, 1),
        A1_time         = np.array([0]), # Only need a placeholder for xp
        boxsize         = 2 * 5.0 * 1e6,
        h_max           = 0.5 * 1e6,
        soft            = 0.09 * 1e6,
        snap_id_type    = "sequence",
    )

    A1900_n65       = Simulation(
        name            = "A1900_n65",
        ic              = init_A1900_n65,
        po              = plot_opt_phodei_init,
        dir_proj        = ut.dir_gihr + "swift/phodei/",
        category        = "phodei_init",
        file_to_SI      = Conversions(1e24, 1e6, 1),
        A1_time         = np.arange(0, 1000 + 1, 100),
        boxsize         = 2 * 500e3,
        h_max           = 0.1 * 1e6,
        soft            = 0.001 * 1e6,
        snap_id_type    = "sequence",
    )

    # Masses
    A1800_n60       = Simulation(
        name            = "A1800_n60",
        copy            = "A1900_n65",
        ic              = init_A1800_n60,
        soft            = 0.0007 * 1e6,
    )
    A1850_n60       = Simulation(
        name            = "A1850_n60",
        copy            = "A1900_n65",
        ic              = init_A1850_n60,
        soft            = 0.001 * 1e6,
    )
    A1950_n65       = Simulation(
        name            = "A1950_n65",
        copy            = "A1900_n65",
        ic              = init_A1950_n65,
        soft            = 0.0015 * 1e6,
    )
    A2000_n65       = Simulation(
        name            = "A2000_n65",
        copy            = "A1900_n65",
        ic              = init_A2000_n65,
        soft            = 0.0021 * 1e6,
    )
    A2050_n70       = Simulation(
        name            = "A2050_n70",
        copy            = "A1900_n65",
        ic              = init_A2050_n70,
        soft            = 0.0021 * 1e6,
    )
    A2100_n70       = Simulation(
        name            = "A2100_n70",
        copy            = "A1900_n65",
        ic              = init_A2100_n70,
        soft            = 0.0032 * 1e6,
    )

    # Spins
    A2000_s030_n65  = Simulation(
        name            = "A2000_s030_n65",
        copy            = "A2000_n65",
        ic              = init_A2000_s030_n65,
    )
    A2000_s036_n65  = Simulation(
        name            = "A2000_s036_n65",
        copy            = "A2000_n65",
        ic              = init_A2000_s036_n65,
    )
    A2000_s047_n65  = Simulation(
        name            = "A2000_s047_n65",
        copy            = "A2000_n65",
        ic              = init_A2000_s047_n65,
    )
    A2000_s086_n65  = Simulation(
        name            = "A2000_s086_n65",
        copy            = "A2000_n65",
        ic              = init_A2000_s086_n65,
    )
    A2000_s170_n65  = Simulation(
        name            = "A2000_s170_n65",
        copy            = "A2000_n65",
        ic              = init_A2000_s170_n65,
    )

    # Differentiated
    A2000c30_n65    = Simulation(
        name            = "A2000c30_n65",
        copy            = "A2000_n65",
        ic              = init_A2000c30_n65,
    )

    # Resolutions
    A2000_n50       = Simulation(
        name            = "A2000_n50",
        copy            = "A1900_n65",
        ic              = init_A2000_n50,
        soft            = 0.0066 * 1e6,
    )
    A2000_n55       = Simulation(
        name            = "A2000_n55",
        copy            = "A1900_n65",
        ic              = init_A2000_n55,
        soft            = 0.0046 * 1e6,
    )
    A2000_n60       = Simulation(
        name            = "A2000_n60",
        copy            = "A1900_n65",
        ic              = init_A2000_n60,
        soft            = 0.0032 * 1e6,
    )
    A2000_n70       = Simulation(
        name            = "A2000_n70",
        copy            = "A1900_n65",
        ic              = init_A2000_n70,
        soft            = 0.0015 * 1e6,
    )

# Encounter simulations
if True:
    # Snapshot times, extra snipshots for animations (#a)
    A1_time_snap = np.arange(1000, 90000 + 1, 1000)
    A1_time_snip = np.arange(100, 90000 + 1, 100)
    A1_time_snap_12h = np.arange(1000, 43300 + 1, 1000)
    A1_time_snip_12h = np.arange(100, 43300 + 1, 100)

    # Periapses and speeds
    Ma_xp_A2000_n65_r11_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            ic_1            = init_Ma_n50,
            ic_2            = init_A2000_n65,
            Fp_init_cond_1  = (Ma_n50, -1),
            Fp_init_cond_2  = (A2000_n65, -1),
            A1_preset_1     = ["xp"],
            A1_preset_2     = ["v=0"],
            q               = 1.1 * Ma.R,
            v_inf           = 0.0e3,
            t_q             = 1.5 * hour_to_s,
            is_centre_mass  = False,
            is_centre_mom   = False,
            m_xp            = Ma.M,
            A1_pos_xp       = [0, 0, 0],
        ),
        po              = plot_opt_phodei,
        dir_proj        = ut.dir_gihr + "swift/phodei/",
        category        = "phodei",
        file_to_SI      = Conversions(1e24, 1e6, 1),
        A1_time         = A1_time_snap,
        boxsize         = 2 * 50 * Ma.R,
        h_max           = 0.042 * 1e6,
        soft            = 0.0021 * 1e6,
        link_len        = "0.0040",
        # A1_link_len     = ["0.0030", "0.0050"],
    )
    Ma_xp_A2000_n65_r12_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_n65_r13_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r13_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.3 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r14_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n65_r14_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.4 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_n65_r15_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r15_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.5 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n65_r16_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.6 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_n65_r17_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r17_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.7 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r18_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n65_r18_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.8 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_n65_r19_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r19_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.9 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r20_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r20_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 2.0 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r22_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r22_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 2.2 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r24_v00     = Simulation(
        name            = "Ma_xp_A2000_n65_r24_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 2.4 * Ma.R,
        ),
    )

    Ma_xp_A2000_n65_r12_v02     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v02",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_n65_r14_v02     = Simulation(
        name            = "Ma_xp_A2000_n65_r14_v02",
        copy            = "Ma_xp_A2000_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v02",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v02     = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v02",
        copy            = "Ma_xp_A2000_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v02",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r18_v02     = Simulation(
        name            = "Ma_xp_A2000_n65_r18_v02",
        copy            = "Ma_xp_A2000_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v02",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r20_v02     = Simulation(
        name            = "Ma_xp_A2000_n65_r20_v02",
        copy            = "Ma_xp_A2000_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v02",
            q               = 2.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_n65_r12_v04     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v04",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_n65_r14_v04     = Simulation(
        name            = "Ma_xp_A2000_n65_r14_v04",
        copy            = "Ma_xp_A2000_n65_r12_v04",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v04",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v04     = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v04",
        copy            = "Ma_xp_A2000_n65_r12_v04",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v04",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r18_v04     = Simulation(
        name            = "Ma_xp_A2000_n65_r18_v04",
        copy            = "Ma_xp_A2000_n65_r12_v04",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v04",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r20_v04     = Simulation(
        name            = "Ma_xp_A2000_n65_r20_v04",
        copy            = "Ma_xp_A2000_n65_r12_v04",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v04",
            q               = 2.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_n65_r12_v06     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v06",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_n65_r14_v06     = Simulation(
        name            = "Ma_xp_A2000_n65_r14_v06",
        copy            = "Ma_xp_A2000_n65_r12_v06",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v06",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v06     = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v06",
        copy            = "Ma_xp_A2000_n65_r12_v06",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v06",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r18_v06     = Simulation(
        name            = "Ma_xp_A2000_n65_r18_v06",
        copy            = "Ma_xp_A2000_n65_r12_v06",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v06",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r20_v06     = Simulation(
        name            = "Ma_xp_A2000_n65_r20_v06",
        copy            = "Ma_xp_A2000_n65_r12_v06",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v06",
            q               = 2.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_n65_r12_v08     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v08",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_n65_r12_v10     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v10",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 1.0e3,
        ),
    )
    Ma_xp_A2000_n65_r12_v12     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v12",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 1.2e3,
        ),
    )
    Ma_xp_A2000_n65_r12_v14     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v14",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 1.4e3,
        ),
    )
    Ma_xp_A2000_n65_r12_v16     = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v16",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 1.6e3,
        ),
    )

    Ma_xp_A2000_n65_r14_v08     = Simulation(
        name            = "Ma_xp_A2000_n65_r14_v08",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.4 * Ma.R,
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_n65_r14_v10     = Simulation(
        name            = "Ma_xp_A2000_n65_r14_v10",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.4 * Ma.R,
            v_inf           = 1.0e3,
        ),
    )

    Ma_xp_A2000_n65_r16_v08     = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v08",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.6 * Ma.R,
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_n65_r16_v10     = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v10",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.6 * Ma.R,
            v_inf           = 1.0e3,
        ),
    )

    Ma_xp_A2000_n65_r18_v08     = Simulation(
        name            = "Ma_xp_A2000_n65_r18_v08",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.8 * Ma.R,
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_n65_r18_v10     = Simulation(
        name            = "Ma_xp_A2000_n65_r18_v10",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            q               = 1.8 * Ma.R,
            v_inf           = 1.0e3,
        ),
    )

    # Spins
    #     Periods and directions
    Ma_xp_A2000_s030z_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s030z_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_s030_n65,
            Fp_init_cond_2  = (A2000_s030_n65, -1),
            A1_preset_2     = ["spin=z"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s036z_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s036z_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_s036_n65,
            Fp_init_cond_2  = (A2000_s036_n65, -1),
            A1_preset_2     = ["spin=z"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s047z_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s047z_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_s047_n65,
            Fp_init_cond_2  = (A2000_s047_n65, -1),
            A1_preset_2     = ["spin=z"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s086z_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s086z_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_s086_n65,
            Fp_init_cond_2  = (A2000_s086_n65, -1),
            A1_preset_2     = ["spin=z"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s170z_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s170z_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_s170_n65,
            Fp_init_cond_2  = (A2000_s170_n65, -1),
            A1_preset_2     = ["spin=z"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )

    Ma_xp_A2000_s030mz_n65_r12_v00  = Simulation(
        name            = "Ma_xp_A2000_s030mz_n65_r12_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s036mz_n65_r12_v00  = Simulation(
        name            = "Ma_xp_A2000_s036mz_n65_r12_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s047mz_n65_r12_v00  = Simulation(
        name            = "Ma_xp_A2000_s047mz_n65_r12_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s086mz_n65_r12_v00  = Simulation(
        name            = "Ma_xp_A2000_s086mz_n65_r12_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s170mz_n65_r12_v00  = Simulation(
        name            = "Ma_xp_A2000_s170mz_n65_r12_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )

    Ma_xp_A2000_s030x_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s030x_n65_r12_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            A1_preset_2     = ["spin=x"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s036x_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s036x_n65_r12_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
            A1_preset_2     = ["spin=x"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s047x_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s047x_n65_r12_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            A1_preset_2     = ["spin=x"],
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000_s086x_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s086x_n65_r12_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s170x_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s170x_n65_r12_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )

    Ma_xp_A2000_s030y_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s030y_n65_r12_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s036y_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s036y_n65_r12_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s047y_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s047y_n65_r12_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s086y_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s086y_n65_r12_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s170y_n65_r12_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s170y_n65_r12_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )

    Ma_xp_A2000_s030z_n65_r16_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s030z_n65_r16_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s036z_n65_r16_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s036z_n65_r16_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s047z_n65_r16_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s086z_n65_r16_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s086z_n65_r16_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s170z_n65_r16_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s170z_n65_r16_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )

    Ma_xp_A2000_s047mz_n65_r16_v00  = Simulation(
        name            = "Ma_xp_A2000_s047mz_n65_r16_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s086mz_n65_r16_v00  = Simulation(
        name            = "Ma_xp_A2000_s086mz_n65_r16_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )
    Ma_xp_A2000_s170mz_n65_r16_v00  = Simulation(
        name            = "Ma_xp_A2000_s170mz_n65_r16_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
            A1_preset_2     = ["spin=-z"],
        ),
    )

    Ma_xp_A2000_s030x_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s030x_n65_r16_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s036x_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s036x_n65_r16_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s047x_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s047x_n65_r16_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s086x_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s086x_n65_r16_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s170x_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s170x_n65_r16_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )

    Ma_xp_A2000_s030y_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s030y_n65_r16_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s036y_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s036y_n65_r16_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s047y_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s047y_n65_r16_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s086y_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s086y_n65_r16_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s170y_n65_r16_v00   = Simulation(
        name            = "Ma_xp_A2000_s170y_n65_r16_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )

    Ma_xp_A2000_s030z_n65_r20_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s030z_n65_r20_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
    )
    Ma_xp_A2000_s036z_n65_r20_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s036z_n65_r20_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s047z_n65_r20_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
    )
    Ma_xp_A2000_s086z_n65_r20_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s086z_n65_r20_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
    )
    Ma_xp_A2000_s170z_n65_r20_v00   = Simulation(#a
        name            = "Ma_xp_A2000_s170z_n65_r20_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_s030x_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s030x_n65_r20_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s036x_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s036x_n65_r20_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r20_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s047x_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s047x_n65_r20_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s086x_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s086x_n65_r20_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r20_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )
    Ma_xp_A2000_s170x_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s170x_n65_r20_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r20_v00",
            A1_preset_2     = ["spin=x"],
        ),
    )

    Ma_xp_A2000_s030y_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s030y_n65_r20_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s036y_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s036y_n65_r20_v00",
        copy            = "Ma_xp_A2000_s036z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r20_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s047y_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s047y_n65_r20_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s086y_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s086y_n65_r20_v00",
        copy            = "Ma_xp_A2000_s086z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r20_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )
    Ma_xp_A2000_s170y_n65_r20_v00   = Simulation(
        name            = "Ma_xp_A2000_s170y_n65_r20_v00",
        copy            = "Ma_xp_A2000_s170z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r20_v00",
            A1_preset_2     = ["spin=y"],
        ),
    )

    Ma_xp_A2000_s030z_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r16_v06",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s036z_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s036z_n65_r16_v06",
        copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036z_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r16_v06",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s086z_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s086z_n65_r16_v06",
        copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086z_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s170z_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s170z_n65_r16_v06",
        copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170z_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )

    Ma_xp_A2000_s047mz_n65_r16_v06  = Simulation(
        name            = "Ma_xp_A2000_s047mz_n65_r16_v06",
        copy            = "Ma_xp_A2000_s047mz_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047mz_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s086mz_n65_r16_v06  = Simulation(
        name            = "Ma_xp_A2000_s086mz_n65_r16_v06",
        copy            = "Ma_xp_A2000_s086mz_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086mz_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s170mz_n65_r16_v06  = Simulation(
        name            = "Ma_xp_A2000_s170mz_n65_r16_v06",
        copy            = "Ma_xp_A2000_s170mz_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170mz_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )

    Ma_xp_A2000_s030x_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s030x_n65_r16_v06",
        copy            = "Ma_xp_A2000_s030x_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030x_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s036x_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s036x_n65_r16_v06",
        copy            = "Ma_xp_A2000_s036x_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036x_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s047x_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s047x_n65_r16_v06",
        copy            = "Ma_xp_A2000_s047x_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047x_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s086x_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s086x_n65_r16_v06",
        copy            = "Ma_xp_A2000_s086x_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086x_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s170x_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s170x_n65_r16_v06",
        copy            = "Ma_xp_A2000_s170x_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170x_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )

    Ma_xp_A2000_s030y_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s030y_n65_r16_v06",
        copy            = "Ma_xp_A2000_s030y_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030y_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s036y_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s036y_n65_r16_v06",
        copy            = "Ma_xp_A2000_s036y_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s036y_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s047y_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s047y_n65_r16_v06",
        copy            = "Ma_xp_A2000_s047y_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047y_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s086y_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s086y_n65_r16_v06",
        copy            = "Ma_xp_A2000_s086y_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s086y_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s170y_n65_r16_v06   = Simulation(
        name            = "Ma_xp_A2000_s170y_n65_r16_v06",
        copy            = "Ma_xp_A2000_s170y_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s170y_n65_r16_v00",
            v_inf           = 0.6e3,
        ),
    )

    #     Periapses and speeds (L_z ~[1, 1/2] L_max)
    Ma_xp_A2000_s030z_n65_r11_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r11_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.1 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v00
    Ma_xp_A2000_s030z_n65_r13_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r13_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.3 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r14_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r14_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r15_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r15_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.5 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r16_v00
    Ma_xp_A2000_s030z_n65_r17_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r17_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.7 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r18_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r18_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r19_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r19_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 1.9 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r20_v00
    Ma_xp_A2000_s030z_n65_r22_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r22_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 2.2 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r24_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r24_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 2.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r26_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r26_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 2.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r28_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r28_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 2.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_s030z_n65_r30_v00   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r30_v00",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            q               = 3.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_s030z_n65_r12_v02   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v02",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v04   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v04",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v06   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v06",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v08   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v08",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v10   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v10",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 1.0e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v12   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v12",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 1.2e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v14   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v14",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 1.4e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r12_v16   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r12_v16",
        copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r12_v00",
            v_inf           = 1.6e3,
        ),
    )

    Ma_xp_A2000_s030z_n65_r16_v02   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r16_v02",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r16_v04   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r16_v04",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r16_v06
    Ma_xp_A2000_s030z_n65_r16_v08   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r16_v08",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r16_v10   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r16_v10",
        copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r16_v00",
            v_inf           = 1.0e3,
        ),
    )

    Ma_xp_A2000_s030z_n65_r20_v02   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r20_v02",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r20_v04   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r20_v04",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r20_v06   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r20_v06",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r20_v08   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r20_v08",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s030z_n65_r20_v10   = Simulation(
        name            = "Ma_xp_A2000_s030z_n65_r20_v10",
        copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s030z_n65_r20_v00",
            v_inf           = 1.0e3,
        ),
    )

    Ma_xp_A2000_s047z_n65_r11_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r11_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.1 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v00
    Ma_xp_A2000_s047z_n65_r13_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r13_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.3 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r14_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r14_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r15_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r15_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.5 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v00
    Ma_xp_A2000_s047z_n65_r17_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r17_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.7 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r18_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r18_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r19_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r19_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 1.9 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v00
    Ma_xp_A2000_s047z_n65_r22_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r22_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 2.2 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r24_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r24_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 2.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r26_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r26_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 2.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r28_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r28_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 2.8 * Ma.R,
        ),
    )
    Ma_xp_A2000_s047z_n65_r30_v00   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r30_v00",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            q               = 3.0 * Ma.R,
        ),
    )

    Ma_xp_A2000_s047z_n65_r12_v02   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v02",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v04   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v04",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v06   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v06",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v08   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v08",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v10   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v10",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 1.0e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v12   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v12",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 1.2e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v14   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v14",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 1.4e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r12_v16   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r12_v16",
        copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r12_v00",
            v_inf           = 1.6e3,
        ),
    )

    Ma_xp_A2000_s047z_n65_r16_v02   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r16_v02",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v04   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r16_v04",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v06
    Ma_xp_A2000_s047z_n65_r16_v08   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r16_v08",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r16_v10   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r16_v10",
        copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r16_v00",
            v_inf           = 1.0e3,
        ),
    )

    Ma_xp_A2000_s047z_n65_r20_v02   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r20_v02",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v04   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r20_v04",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v06   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r20_v06",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v08   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r20_v08",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000_s047z_n65_r20_v10   = Simulation(
        name            = "Ma_xp_A2000_s047z_n65_r20_v10",
        copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_s047z_n65_r20_v00",
            v_inf           = 1.0e3,
        ),
    )

    # Masses
    Ma_xp_A1800_n60_r12_v00     = Simulation(
        name            = "Ma_xp_A1800_n60_r12_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            ic_2            = init_A1800_n60,
            Fp_init_cond_2  = (A1800_n60, -1),
            q               = 1.2 * Ma.R,
            v_inf           = 0.0e3,
        ),
        h_max           = 0.013 * 1e6,
        soft            = 0.0007 * 1e6,
        link_len        = "0.0013",
        # A1_link_len     = ["0.0016"], #["0.0001", "0.0016"],
    )
    Ma_xp_A1850_n60_r12_v00     = Simulation(
        name            = "Ma_xp_A1850_n60_r12_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            ic_2            = init_A1850_n60,
            Fp_init_cond_2  = (A1850_n60, -1),
            q               = 1.2 * Ma.R,
            v_inf           = 0.0e3,
        ),
        h_max           = 0.019 * 1e6,
        soft            = 0.001 * 1e6,
        link_len        = "0.0019",
        # A1_link_len     = ["0.0014", "0.0024"],
    )
    Ma_xp_A1900_n65_r12_v00     = Simulation(
        name            = "Ma_xp_A1900_n65_r11_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r11_v00",
            ic_2            = init_A1900_n65,
            Fp_init_cond_2  = (A1900_n65, -1),
            q               = 1.2 * Ma.R,
            v_inf           = 0.0e3,
        ),
        h_max           = 0.05 * 1e6,
        soft            = 0.001 * 1e6,
        link_len        = "0.0019",
        # A1_link_len     = ["0.0014", "0.0024"],
    )
    Ma_xp_A1950_n65_r12_v00     = Simulation(
        name            = "Ma_xp_A1950_n65_r12_v00",
        copy            = "Ma_xp_A1800_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1800_n60_r12_v00",
            ic_2            = init_A1950_n65,
            Fp_init_cond_2  = (A1950_n65, -1),
        ),
        h_max           = 0.029 * 1e6,
        soft            = 0.0015 * 1e6,
        link_len        = "0.0027",
        # A1_link_len     = ["0.0020", "0.0034"],
    )
    Ma_xp_A2000_n65_r12_v00
    Ma_xp_A2050_n70_r12_v00     = Simulation(
        name            = "Ma_xp_A2050_n70_r12_v00",
        copy            = "Ma_xp_A1800_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1800_n60_r12_v00",
            ic_2            = init_A2050_n70,
            Fp_init_cond_2  = (A2050_n70, -1),
        ),
        h_max           = 0.042 * 1e6,
        soft            = 0.0021 * 1e6,
        link_len        = "0.0040",
        # A1_link_len     = ["0.0030", "0.0050"],
    )
    Ma_xp_A2100_n70_r12_v00     = Simulation(
        name            = "Ma_xp_A2100_n70_r12_v00",
        copy            = "Ma_xp_A1800_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1800_n60_r12_v00",
            ic_2            = init_A2100_n70,
            Fp_init_cond_2  = (A2100_n70, -1),
        ),
        h_max           = 0.062 * 1e6,
        soft            = 0.0032 * 1e6,
        link_len        = "0.0059",
        # A1_link_len     = ["0.0044", "0.0074"],
    )

    Ma_xp_A1800_n60_r16_v00     = Simulation(
        name            = "Ma_xp_A1800_n60_r16_v00",
        copy            = "Ma_xp_A1800_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1800_n60_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A1850_n60_r16_v00     = Simulation(
        name            = "Ma_xp_A1850_n60_r16_v00",
        copy            = "Ma_xp_A1850_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1850_n60_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A1900_n65_r16_v00     = Simulation(
        name            = "Ma_xp_A1900_n65_r16_v00",
        copy            = "Ma_xp_A1900_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1900_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A1950_n65_r16_v00     = Simulation(
        name            = "Ma_xp_A1950_n65_r16_v00",
        copy            = "Ma_xp_A1950_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A1950_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v00
    Ma_xp_A2050_n70_r16_v00     = Simulation(
        name            = "Ma_xp_A2050_n70_r16_v00",
        copy            = "Ma_xp_A2050_n70_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2050_n70_r12_v00",
            q               = 1.6 * Ma.R,
        )
    )
    Ma_xp_A2100_n70_r16_v00     = Simulation(
        name            = "Ma_xp_A2100_n70_r16_v00",
        copy            = "Ma_xp_A2100_n70_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2100_n70_r12_v00",
            q               = 1.6 * Ma.R,
        )
    )

    # Differentiated
    Ma_xp_A2000c30_n65_r12_v00  = Simulation(#a
        name            = "Ma_xp_A2000c30_n65_r12_v00",
        copy            = "Ma_xp_A2000_n65_r11_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000c30_n65,
            Fp_init_cond_2  = (A2000c30_n65, -1),
            q               = 1.2 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000c30_n65_r14_v00  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r14_v00",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000c30_n65_r16_v00  = Simulation(#a
        name            = "Ma_xp_A2000c30_n65_r16_v00",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 1.6 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000c30_n65_r18_v00  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r18_v00",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 1.8 * Ma.R,
        ),
    )
    Ma_xp_A2000c30_n65_r20_v00  = Simulation(#a
        name            = "Ma_xp_A2000c30_n65_r20_v00",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 2.0 * Ma.R,
        ),
        A1_time_snap    = A1_time_snap,
        A1_time_snip    = A1_time_snip,
    )
    Ma_xp_A2000c30_n65_r22_v00  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r22_v00",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 2.2 * Ma.R,
        ),
    )

    Ma_xp_A2000c30_n65_r12_v02  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r12_v02",
        copy            = "Ma_xp_A2000c30_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v00",
            q               = 1.2 * Ma.R,
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000c30_n65_r12_v04  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r12_v04",
        copy            = "Ma_xp_A2000c30_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v02",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000c30_n65_r12_v06  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r12_v06",
        copy            = "Ma_xp_A2000c30_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v02",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000c30_n65_r12_v08  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r12_v08",
        copy            = "Ma_xp_A2000c30_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v02",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000c30_n65_r12_v10  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r12_v10",
        copy            = "Ma_xp_A2000c30_n65_r12_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r12_v02",
            v_inf           = 01.0e3,
        ),
    )

    Ma_xp_A2000c30_n65_r16_v02  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r16_v02",
        copy            = "Ma_xp_A2000c30_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r16_v00",
            q               = 1.6 * Ma.R,
            v_inf           = 0.2e3,
        ),
    )
    Ma_xp_A2000c30_n65_r16_v04  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r16_v04",
        copy            = "Ma_xp_A2000c30_n65_r16_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r16_v02",
            v_inf           = 0.4e3,
        ),
    )
    Ma_xp_A2000c30_n65_r16_v06  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r16_v06",
        copy            = "Ma_xp_A2000c30_n65_r16_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r16_v02",
            v_inf           = 0.6e3,
        ),
    )
    Ma_xp_A2000c30_n65_r16_v08  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r16_v08",
        copy            = "Ma_xp_A2000c30_n65_r16_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r16_v02",
            v_inf           = 0.8e3,
        ),
    )
    Ma_xp_A2000c30_n65_r16_v10  = Simulation(
        name            = "Ma_xp_A2000c30_n65_r16_v10",
        copy            = "Ma_xp_A2000c30_n65_r16_v02",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000c30_n65_r16_v02",
            v_inf           = 01.0e3,
        ),
    )

    # Resolutions
    Ma_xp_A2000_n50_r12_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n50_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_n50,
            Fp_init_cond_2  = (A2000_n50, -1),
        ),
        h_max           = 0.13 * 1e6,
        soft            = 0.0066 * 1e6,
        link_len        = "0.0126",
    )
    Ma_xp_A2000_n55_r12_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n55_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_n55,
            Fp_init_cond_2  = (A2000_n55, -1),
        ),
        h_max           = 0.091 * 1e6,
        soft            = 0.0046 * 1e6,
        link_len        = "0.0086",
    )
    Ma_xp_A2000_n60_r12_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n60_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_n60,
            Fp_init_cond_2  = (A2000_n60, -1),
        ),
        h_max           = 0.062 * 1e6,
        soft            = 0.0032 * 1e6,
        link_len        = "0.0059",
    )
    Ma_xp_A2000_n65_r12_v00
    Ma_xp_A2000_n70_r12_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n70_r12_v00",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            ic_2            = init_A2000_n70,
            Fp_init_cond_2  = (A2000_n70, -1),
        ),
        dir_proj        = ut.dir_gihr8 + "swift/phodei/",
        h_max           = 0.029 * 1e6,
        soft            = 0.0015 * 1e6,
        link_len        = "0.0020",
    )

    Ma_xp_A2000_n50_r14_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n50_r14_v00",
        copy            = "Ma_xp_A2000_n50_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n50_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n55_r14_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n55_r14_v00",
        copy            = "Ma_xp_A2000_n55_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n55_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n60_r14_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n60_r14_v00",
        copy            = "Ma_xp_A2000_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n60_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r14_v00
    Ma_xp_A2000_n70_r14_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n70_r14_v00",
        copy            = "Ma_xp_A2000_n70_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n70_r12_v00",
            q               = 1.4 * Ma.R,
        ),
    )

    Ma_xp_A2000_n50_r16_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n50_r16_v00",
        copy            = "Ma_xp_A2000_n50_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n50_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n55_r16_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n55_r16_v00",
        copy            = "Ma_xp_A2000_n55_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n55_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n60_r16_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n60_r16_v00",
        copy            = "Ma_xp_A2000_n60_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n60_r12_v00",
            q               = 1.6 * Ma.R,
        ),
    )
    Ma_xp_A2000_n65_r16_v00
    Ma_xp_A2000_n70_r16_v00     = Simulation(#a
        name            = "Ma_xp_A2000_n70_r16_v00",
        copy            = "Ma_xp_A2000_n70_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n70_r12_v00",
            q               = 1.6 * Ma.R,
        ),
        A1_time_snap    = np.arange(1000, 90000 + 1, 1000),
        A1_time_snip    = np.arange(100, 90000 + 1, 100),
        Di_po_edit      = {
            "func_snapshot_add_misc"    : ut_ph.add_panel_text_disr_snaps,
        }
    )

    # Reoriented repeats
    Ma_xp_A2000_n65_r12_v00_r1  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r1",
        copy            = "Ma_xp_A2000_n65_r12_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00",
            A1_preset_2     = ["v=0", "rotate=z30"],
        ),
        A1_time         = A1_time_snap,
    )
    Ma_xp_A2000_n65_r12_v00_r2  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r2",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z60"],
        ),
    )
    Ma_xp_A2000_n65_r12_v00_r3  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r3",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z90"],
        ),
    )
    Ma_xp_A2000_n65_r12_v00_r4  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r4",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z120"],
        ),
    )
    Ma_xp_A2000_n65_r12_v00_r5  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r5",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z150"],
        ),
    )
    Ma_xp_A2000_n65_r12_v00_r6  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r6",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z180"],
        ),
    )
    Ma_xp_A2000_n65_r12_v00_r7  = Simulation(
        name            = "Ma_xp_A2000_n65_r12_v00_r7",
        copy            = "Ma_xp_A2000_n65_r12_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r12_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z210"],
        ),
    )

    Ma_xp_A2000_n65_r16_v00_r1  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r1",
        copy            = "Ma_xp_A2000_n65_r16_v00",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00",
            A1_preset_2     = ["v=0", "rotate=z30"],
        ),
        A1_time         = A1_time_snap,
    )
    Ma_xp_A2000_n65_r16_v00_r2  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r2",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z60"],
        ),
    )
    Ma_xp_A2000_n65_r16_v00_r3  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r3",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z90"],
        ),
    )
    Ma_xp_A2000_n65_r16_v00_r4  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r4",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z120"],
        ),
    )
    Ma_xp_A2000_n65_r16_v00_r5  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r5",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z150"],
        ),
    )
    Ma_xp_A2000_n65_r16_v00_r6  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r6",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z180"],
        ),
    )
    Ma_xp_A2000_n65_r16_v00_r7  = Simulation(
        name            = "Ma_xp_A2000_n65_r16_v00_r7",
        copy            = "Ma_xp_A2000_n65_r16_v00_r1",
        impact          = ImpactInitCond(
            copy            = "Ma_xp_A2000_n65_r16_v00_r1",
            A1_preset_2     = ["v=0", "rotate=z210"],
        ),
    )

# Simulation sets
if True:
    # Periapses and speeds
    set_Ma_xp_A2000_n65_r___v00     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v00",
        A1_sim          = [
            Ma_xp_A2000_n65_r11_v00,
            Ma_xp_A2000_n65_r12_v00,
            Ma_xp_A2000_n65_r13_v00,
            Ma_xp_A2000_n65_r14_v00,
            Ma_xp_A2000_n65_r15_v00,
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2000_n65_r17_v00,
            Ma_xp_A2000_n65_r18_v00,
            Ma_xp_A2000_n65_r19_v00,
            Ma_xp_A2000_n65_r20_v00,
            Ma_xp_A2000_n65_r22_v00,
            Ma_xp_A2000_n65_r24_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_n65_r___v02     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v02",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v02,
            Ma_xp_A2000_n65_r14_v02,
            Ma_xp_A2000_n65_r16_v02,
            Ma_xp_A2000_n65_r18_v02,
            Ma_xp_A2000_n65_r20_v02,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_n65_r___v04     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v04",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v04,
            Ma_xp_A2000_n65_r14_v04,
            Ma_xp_A2000_n65_r16_v04,
            Ma_xp_A2000_n65_r18_v04,
            Ma_xp_A2000_n65_r20_v04,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_n65_r___v06     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v06",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v06,
            Ma_xp_A2000_n65_r14_v06,
            Ma_xp_A2000_n65_r16_v06,
            Ma_xp_A2000_n65_r18_v06,
            Ma_xp_A2000_n65_r20_v06,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_n65_r___v08     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v08",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v08,
            Ma_xp_A2000_n65_r14_v08,
            Ma_xp_A2000_n65_r16_v08,
            Ma_xp_A2000_n65_r18_v08,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_n65_r___v10     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v08",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v10,
            Ma_xp_A2000_n65_r14_v10,
            Ma_xp_A2000_n65_r16_v10,
            Ma_xp_A2000_n65_r18_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )

    set_Ma_xp_A2000_n65_r12_v___1   = SimSet(
        name            = "set_Ma_xp_A2000_n65_r12_v___1",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v02,
            Ma_xp_A2000_n65_r12_v04,
            Ma_xp_A2000_n65_r12_v06,
            Ma_xp_A2000_n65_r12_v08,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r12_v___2   = SimSet(
        name            = "set_Ma_xp_A2000_n65_r12_v___2",
        A1_sim          = [
            Ma_xp_A2000_n65_r12_v10,
            Ma_xp_A2000_n65_r12_v12,
            Ma_xp_A2000_n65_r12_v14,
            Ma_xp_A2000_n65_r12_v16,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r12_v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r12_v__",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_n65_r12_v___1.A1_sim,
            set_Ma_xp_A2000_n65_r12_v___2.A1_sim,
        )),
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r14_v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r14_v__",
        A1_sim          = [
            Ma_xp_A2000_n65_r14_v00,
            Ma_xp_A2000_n65_r14_v02,
            Ma_xp_A2000_n65_r14_v04,
            Ma_xp_A2000_n65_r14_v06,
            Ma_xp_A2000_n65_r14_v08,
            Ma_xp_A2000_n65_r14_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r16_v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r16_v__",
        A1_sim          = [
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v02,
            Ma_xp_A2000_n65_r16_v04,
            Ma_xp_A2000_n65_r16_v06,
            Ma_xp_A2000_n65_r16_v08,
            Ma_xp_A2000_n65_r16_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r18_v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r18_v__",
        A1_sim          = [
            Ma_xp_A2000_n65_r18_v00,
            Ma_xp_A2000_n65_r18_v02,
            Ma_xp_A2000_n65_r18_v04,
            Ma_xp_A2000_n65_r18_v06,
            Ma_xp_A2000_n65_r18_v08,
            Ma_xp_A2000_n65_r18_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_n65_r20_v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r20_v__",
        A1_sim          = [
            Ma_xp_A2000_n65_r20_v00,
            Ma_xp_A2000_n65_r20_v02,
            Ma_xp_A2000_n65_r20_v04,
            Ma_xp_A2000_n65_r20_v06,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )

    set_Ma_xp_A2000_n65_r___v0_     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v0_",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_n65_r___v02.A1_sim,
            set_Ma_xp_A2000_n65_r___v04.A1_sim,
            set_Ma_xp_A2000_n65_r___v06.A1_sim,
            set_Ma_xp_A2000_n65_r___v08.A1_sim,
            set_Ma_xp_A2000_n65_r___v10.A1_sim,
        )),
    )
    set_Ma_xp_A2000_n65_r___v__     = SimSet(
        name            = "set_Ma_xp_A2000_n65_r___v__",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_n65_r___v0_.A1_sim,
            set_Ma_xp_A2000_n65_r12_v__.A1_sim,
            set_Ma_xp_A2000_n65_r14_v__.A1_sim,
            set_Ma_xp_A2000_n65_r16_v__.A1_sim,
            set_Ma_xp_A2000_n65_r18_v__.A1_sim,
            set_Ma_xp_A2000_n65_r20_v__.A1_sim,
        )),
        A1_colour       = [cmap_rbow],
    )

    # Spins
    #     Periods and directions
    set_Ma_xp_A2000_s___z_n65_r12_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r12_v00,
            Ma_xp_A2000_s036z_n65_r12_v00,
            Ma_xp_A2000_s047z_n65_r12_v00,
            Ma_xp_A2000_s086z_n65_r12_v00,
            Ma_xp_A2000_s170z_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___mz_n65_r12_v00  = SimSet(
        name            = "set_Ma_xp_A2000_s___mz_n65_r12_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030mz_n65_r12_v00,
            Ma_xp_A2000_s036mz_n65_r12_v00,
            Ma_xp_A2000_s047mz_n65_r12_v00,
            Ma_xp_A2000_s086mz_n65_r12_v00,
            Ma_xp_A2000_s170mz_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
        ],
        legend_title    = r"$L_{-z}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___x_n65_r12_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___x_n65_r12_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030x_n65_r12_v00,
            Ma_xp_A2000_s036x_n65_r12_v00,
            Ma_xp_A2000_s047x_n65_r12_v00,
            Ma_xp_A2000_s086x_n65_r12_v00,
            Ma_xp_A2000_s170x_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
        ],
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___y_n65_r12_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___y_n65_r12_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030y_n65_r12_v00,
            Ma_xp_A2000_s036y_n65_r12_v00,
            Ma_xp_A2000_s047y_n65_r12_v00,
            Ma_xp_A2000_s086y_n65_r12_v00,
            Ma_xp_A2000_s170y_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
        ],
        legend_title    = r"$L_{y}$ $(L_{\rm max})$",
    )

    set_Ma_xp_A2000_s___z_n65_r16_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___z_n65_r16_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r16_v00,
            Ma_xp_A2000_s036z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s086z_n65_r16_v00,
            Ma_xp_A2000_s170z_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
        ],
    )
    set_Ma_xp_A2000_s___mz_n65_r16_v00  = SimSet(
        name            = "set_Ma_xp_A2000_s___mz_n65_r16_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v00",
        A1_sim          = [
            Ma_xp_A2000_s047mz_n65_r16_v00,
            Ma_xp_A2000_s086mz_n65_r16_v00,
            Ma_xp_A2000_s170mz_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
        ],
        A1_colour       = [cmap_rbow],
        legend_title    = r"$L_{-z}$ $(L_{\rm max})$",
        A1_label        = ["1/2", "1/4", "1/8", "0"],  # L_max
    )
    set_Ma_xp_A2000_s___x_n65_r16_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___x_n65_r16_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v00",
        A1_sim          = [
            Ma_xp_A2000_s030x_n65_r16_v00,
            Ma_xp_A2000_s036x_n65_r16_v00,
            Ma_xp_A2000_s047x_n65_r16_v00,
            Ma_xp_A2000_s086x_n65_r16_v00,
            Ma_xp_A2000_s170x_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
        ],
        A1_colour       = [cmap_rbow_r],
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___y_n65_r16_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___y_n65_r16_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v00",
        A1_sim          = [
            Ma_xp_A2000_s030y_n65_r16_v00,
            Ma_xp_A2000_s036y_n65_r16_v00,
            Ma_xp_A2000_s047y_n65_r16_v00,
            Ma_xp_A2000_s086y_n65_r16_v00,
            Ma_xp_A2000_s170y_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
        ],
        A1_colour       = [cmap_rbow_r],
        legend_title    = r"$L_{y}$ $(L_{\rm max})$",
    )

    set_Ma_xp_A2000_s___z_n65_r20_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___z_n65_r20_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r20_v00,
            Ma_xp_A2000_s036z_n65_r20_v00,
            Ma_xp_A2000_s047z_n65_r20_v00,
            Ma_xp_A2000_s086z_n65_r20_v00,
            Ma_xp_A2000_s170z_n65_r20_v00,
            Ma_xp_A2000_n65_r20_v00,
        ],
    )
    set_Ma_xp_A2000_s___x_n65_r20_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___x_n65_r20_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r20_v00",
        A1_sim          = [
            Ma_xp_A2000_s030x_n65_r20_v00,
            Ma_xp_A2000_s036x_n65_r20_v00,
            Ma_xp_A2000_s047x_n65_r20_v00,
            Ma_xp_A2000_s086x_n65_r20_v00,
            Ma_xp_A2000_s170x_n65_r20_v00,
            Ma_xp_A2000_n65_r20_v00,
        ],
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___y_n65_r20_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s___y_n65_r20_v00",
        copy            = "set_Ma_xp_A2000_s___z_n65_r20_v00",
        A1_sim          = [
            Ma_xp_A2000_s030y_n65_r20_v00,
            Ma_xp_A2000_s036y_n65_r20_v00,
            Ma_xp_A2000_s047y_n65_r20_v00,
            Ma_xp_A2000_s086y_n65_r20_v00,
            Ma_xp_A2000_s170y_n65_r20_v00,
            Ma_xp_A2000_n65_r20_v00,
        ],
        legend_title    = r"$L_{y}$ $(L_{\rm max})$",
    )

    set_Ma_xp_A2000_s___z_n65_r16_v06   = SimSet(
        name            = "set_Ma_xp_A2000_s___z_n65_r16_v06",
        copy            = "set_Ma_xp_A2000_s___z_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r16_v06,
            Ma_xp_A2000_s036z_n65_r16_v06,
            Ma_xp_A2000_s047z_n65_r16_v06,
            Ma_xp_A2000_s086z_n65_r16_v06,
            Ma_xp_A2000_s170z_n65_r16_v06,
            Ma_xp_A2000_n65_r16_v06,
        ],
    )
    set_Ma_xp_A2000_s___mz_n65_r16_v06  = SimSet(
        name            = "set_Ma_xp_A2000_s___mz_n65_r16_v06",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v06",
        A1_sim          = [
            Ma_xp_A2000_s047mz_n65_r16_v06,
            Ma_xp_A2000_s086mz_n65_r16_v06,
            Ma_xp_A2000_s170mz_n65_r16_v06,
            Ma_xp_A2000_n65_r16_v06,
        ],
        A1_colour       = [cmap_rbow],
        legend_title    = r"$L_{-z}$ $(L_{\rm max})$",
        A1_label        = ["1/2", "1/4", "1/8", "0"],  # L_max
    )
    set_Ma_xp_A2000_s___x_n65_r16_v06   = SimSet(
        name            = "set_Ma_xp_A2000_s___x_n65_r16_v06",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v06",
        A1_sim          = [
            Ma_xp_A2000_s030x_n65_r16_v06,
            Ma_xp_A2000_s036x_n65_r16_v06,
            Ma_xp_A2000_s047x_n65_r16_v06,
            Ma_xp_A2000_s086x_n65_r16_v06,
            Ma_xp_A2000_s170x_n65_r16_v06,
            Ma_xp_A2000_n65_r16_v06,
        ],
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___y_n65_r16_v06   = SimSet(
        name            = "set_Ma_xp_A2000_s___y_n65_r16_v06",
        copy            = "set_Ma_xp_A2000_s___z_n65_r16_v06",
        A1_sim          = [
            Ma_xp_A2000_s030y_n65_r16_v06,
            Ma_xp_A2000_s036y_n65_r16_v06,
            Ma_xp_A2000_s047y_n65_r16_v06,
            Ma_xp_A2000_s086y_n65_r16_v06,
            Ma_xp_A2000_s170y_n65_r16_v06,
            Ma_xp_A2000_n65_r16_v06,
        ],
        legend_title    = r"$L_{y}$ $(L_{\rm max})$",
    )

    set_Ma_xp_A2000_s_____n65_r12_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s_____n65_r12_v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s___z_n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s___mz_n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s___x_n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s___y_n65_r12_v00.A1_sim,
        )),
    )
    set_Ma_xp_A2000_s_____n65_r16_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s_____n65_r16_v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s___z_n65_r16_v00.A1_sim,
            set_Ma_xp_A2000_s___mz_n65_r16_v00.A1_sim,
            set_Ma_xp_A2000_s___x_n65_r16_v00.A1_sim,
            set_Ma_xp_A2000_s___y_n65_r16_v00.A1_sim,
        )),
    )
    set_Ma_xp_A2000_s_____n65_r20_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s_____n65_r20_v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s___z_n65_r20_v00.A1_sim,
            set_Ma_xp_A2000_s___x_n65_r20_v00.A1_sim,
            set_Ma_xp_A2000_s___y_n65_r20_v00.A1_sim,
        )),
    )
    set_Ma_xp_A2000_s_____n65_r16_v06   = SimSet(
        name            = "set_Ma_xp_A2000_s_____n65_r16_v06",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s___z_n65_r16_v06.A1_sim,
            set_Ma_xp_A2000_s___mz_n65_r16_v06.A1_sim,
            set_Ma_xp_A2000_s___x_n65_r16_v06.A1_sim,
            set_Ma_xp_A2000_s___y_n65_r16_v06.A1_sim,
        )),
    )
    set_Ma_xp_A2000_s_____n65_r___v0_   = SimSet(
        name            = "set_Ma_xp_A2000_s_____n65_r___v0_",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s_____n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s_____n65_r16_v00.A1_sim,
            set_Ma_xp_A2000_s_____n65_r20_v00.A1_sim,
            # set_Ma_xp_A2000_s_____n65_r16_v06.A1_sim,
        )),
    )

    set_Ma_xp_A2000_s___zmz_n65_r12_v00 = SimSet(
        name            = "set_Ma_xp_A2000_s___zmz_n65_r12_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r12_v00,
            Ma_xp_A2000_s036z_n65_r12_v00,
            Ma_xp_A2000_s047z_n65_r12_v00,
            Ma_xp_A2000_s086z_n65_r12_v00,
            Ma_xp_A2000_s170z_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
            Ma_xp_A2000_s170mz_n65_r12_v00,
            Ma_xp_A2000_s086mz_n65_r12_v00,
            Ma_xp_A2000_s047mz_n65_r12_v00,
            Ma_xp_A2000_s036mz_n65_r12_v00,
            Ma_xp_A2000_s030mz_n65_r12_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = [
            "1", "3/4", "1/2", "1/4", "1/8",
            "0",
            r"$-$1/8", r"$-$1/4", r"$-$1/2", r"$-$3/4", r"$-$1",
        ],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
    )
    set_Ma_xp_A2000_s___zmz_n65_r16_v00 = SimSet(
        name            = "set_Ma_xp_A2000_s___zmz_n65_r16_v00",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r16_v00,
            Ma_xp_A2000_s036z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s086z_n65_r16_v00,
            Ma_xp_A2000_s170z_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2000_s170mz_n65_r16_v00,
            Ma_xp_A2000_s086mz_n65_r16_v00,
            Ma_xp_A2000_s047mz_n65_r16_v00,
        ],
        A1_colour       = [cmap_rbow_r],
        A1_label        = [
            "1", "3/4", "1/2", "1/4", "1/8",
            "0",
            r"$-$1/8", r"$-$1/4", r"$-$1/2",
        ],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
    )

    set_Ma_xp_A2000_s047__n65_r16_v00   = SimSet(
        name            = "set_Ma_xp_A2000_s047__n65_r16_v00",
        A1_sim          = [
            Ma_xp_A2000_n70_r16_v00,
            Ma_xp_A2000_s047x_n65_r16_v00,
            Ma_xp_A2000_s047y_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
        ],
        A1_colour       = ["k", A1_c[0], A1_c[1], A1_c[2]],
        A1_label        = [r"$0$", r"$x$", r"$y$", r"$z$"],
        legend_title    = r"Spin Axis",
    )

    set_Ma_xp_A2000_s___zmz_n65_r12r16_v00  = SimSet(
        name            = "set_Ma_xp_A2000_s___zmz_n65_r12r16_v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s___zmz_n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s___zmz_n65_r16_v00.A1_sim,
        )),
    )

    #   Periapses and speeds (~[1, 1/2] L_max)
    set_Ma_xp_A2000_s030z_n65_r___v00_1 = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r___v00_1",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r11_v00,
            Ma_xp_A2000_s030z_n65_r12_v00,
            Ma_xp_A2000_s030z_n65_r13_v00,
            Ma_xp_A2000_s030z_n65_r14_v00,
            Ma_xp_A2000_s030z_n65_r15_v00,
            Ma_xp_A2000_s030z_n65_r16_v00,
            Ma_xp_A2000_s030z_n65_r17_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_s030z_n65_r___v00_2 = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r___v00_2",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r18_v00,
            Ma_xp_A2000_s030z_n65_r19_v00,
            Ma_xp_A2000_s030z_n65_r20_v00,
            Ma_xp_A2000_s030z_n65_r22_v00,
            Ma_xp_A2000_s030z_n65_r24_v00,
            Ma_xp_A2000_s030z_n65_r26_v00,
            Ma_xp_A2000_s030z_n65_r28_v00,
            Ma_xp_A2000_s030z_n65_r30_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_s030z_n65_r12_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r12_v__",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r12_v00,
            Ma_xp_A2000_s030z_n65_r12_v02,
            Ma_xp_A2000_s030z_n65_r12_v04,
            Ma_xp_A2000_s030z_n65_r12_v06,
            Ma_xp_A2000_s030z_n65_r12_v08,
            Ma_xp_A2000_s030z_n65_r12_v10,
            Ma_xp_A2000_s030z_n65_r12_v12,
            Ma_xp_A2000_s030z_n65_r12_v14,
            Ma_xp_A2000_s030z_n65_r12_v16,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_s030z_n65_r16_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r16_v__",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r16_v00,
            Ma_xp_A2000_s030z_n65_r16_v02,
            Ma_xp_A2000_s030z_n65_r16_v04,
            Ma_xp_A2000_s030z_n65_r16_v06,
            Ma_xp_A2000_s030z_n65_r16_v08,
            Ma_xp_A2000_s030z_n65_r16_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_s030z_n65_r20_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r20_v__",
        A1_sim          = [
            Ma_xp_A2000_s030z_n65_r20_v00,
            Ma_xp_A2000_s030z_n65_r20_v02,
            Ma_xp_A2000_s030z_n65_r20_v04,
            Ma_xp_A2000_s030z_n65_r20_v06,
            Ma_xp_A2000_s030z_n65_r20_v08,
            Ma_xp_A2000_s030z_n65_r20_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )

    set_Ma_xp_A2000_s047z_n65_r___v00_1 = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r___v00_1",
        A1_sim          = [
            Ma_xp_A2000_s047z_n65_r11_v00,
            Ma_xp_A2000_s047z_n65_r12_v00,
            Ma_xp_A2000_s047z_n65_r13_v00,
            Ma_xp_A2000_s047z_n65_r14_v00,
            Ma_xp_A2000_s047z_n65_r15_v00,
            Ma_xp_A2000_s047z_n65_r17_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_s047z_n65_r___v00_2 = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r___v00_2",
        A1_sim          = [
            Ma_xp_A2000_s047z_n65_r18_v00,
            Ma_xp_A2000_s047z_n65_r19_v00,
            Ma_xp_A2000_s047z_n65_r20_v00,
            Ma_xp_A2000_s047z_n65_r22_v00,
            Ma_xp_A2000_s047z_n65_r24_v00,
            Ma_xp_A2000_s047z_n65_r26_v00,
            Ma_xp_A2000_s047z_n65_r28_v00,
            Ma_xp_A2000_s047z_n65_r30_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_s047z_n65_r12_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r12_v__",
        A1_sim          = [
            Ma_xp_A2000_s047z_n65_r12_v00,
            Ma_xp_A2000_s047z_n65_r12_v02,
            Ma_xp_A2000_s047z_n65_r12_v04,
            Ma_xp_A2000_s047z_n65_r12_v06,
            Ma_xp_A2000_s047z_n65_r12_v08,
            Ma_xp_A2000_s047z_n65_r12_v10,
            Ma_xp_A2000_s047z_n65_r12_v12,
            Ma_xp_A2000_s047z_n65_r12_v14,
            Ma_xp_A2000_s047z_n65_r12_v16,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_s047z_n65_r16_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r16_v__",
        A1_sim          = [
            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v02,
            Ma_xp_A2000_s047z_n65_r16_v04,
            Ma_xp_A2000_s047z_n65_r16_v06,
            Ma_xp_A2000_s047z_n65_r16_v08,
            Ma_xp_A2000_s047z_n65_r16_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000_s047z_n65_r20_v__   = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r20_v__",
        A1_sim          = [
            Ma_xp_A2000_s047z_n65_r20_v00,
            Ma_xp_A2000_s047z_n65_r20_v02,
            Ma_xp_A2000_s047z_n65_r20_v04,
            Ma_xp_A2000_s047z_n65_r20_v06,
            Ma_xp_A2000_s047z_n65_r20_v08,
            Ma_xp_A2000_s047z_n65_r20_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )

    set_Ma_xp_A2000_s030z_n65_r___v00   = SimSet(
        name            = "set_Ma_xp_A2000_s030z_n65_r___v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s030z_n65_r___v00_1.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r___v00_2.A1_sim,
        )),
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )
    set_Ma_xp_A2000_s047z_n65_r___v00   = SimSet(
        name            = "set_Ma_xp_A2000_s047z_n65_r___v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_s047z_n65_r___v00_1.A1_sim,
            set_Ma_xp_A2000_s047z_n65_r___v00_2.A1_sim,
        )),
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=q"],
    )

    # Masses
    set_Ma_xp_A_____n___r12_v00     = SimSet(
        name            = "set_Ma_xp_A_____n___r12_v00",
        A1_sim          = [
            Ma_xp_A1800_n60_r12_v00,
            Ma_xp_A1850_n60_r12_v00,
            Ma_xp_A1900_n65_r12_v00,
            Ma_xp_A1950_n65_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
            Ma_xp_A2050_n70_r12_v00,
            Ma_xp_A2100_n70_r12_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=M_i"],
    )
    set_Ma_xp_A_____n___r16_v00     = SimSet(
        name            = "set_Ma_xp_A_____n___r16_v00",
        copy            = "set_Ma_xp_A_____n___r12_v00",
        A1_sim          = [
            Ma_xp_A1800_n60_r16_v00,
            Ma_xp_A1850_n60_r16_v00,
            Ma_xp_A1900_n65_r16_v00,
            Ma_xp_A1950_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2050_n70_r16_v00,
            Ma_xp_A2100_n70_r16_v00,
        ],
    )
    set_Ma_xp_A_____n___r12r16_v00  = SimSet(
        name            = "set_Ma_xp_A_____n___r12r16_v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A_____n___r12_v00.A1_sim,
            set_Ma_xp_A_____n___r16_v00.A1_sim,
        )),
    )

    # Differentiated
    set_Ma_xp_A2000c30_n65_r___v00      = SimSet(
        name            = "set_Ma_xp_A2000c30_n65_r___v00",
        A1_sim          = [
            Ma_xp_A2000c30_n65_r12_v00,
            Ma_xp_A2000c30_n65_r14_v00,
            Ma_xp_A2000c30_n65_r16_v00,
            Ma_xp_A2000c30_n65_r18_v00,
            Ma_xp_A2000c30_n65_r20_v00,
            Ma_xp_A2000c30_n65_r22_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=r_p"],
    )
    set_Ma_xp_A2000c30_n65_r12_v__      = SimSet(
        name            = "set_Ma_xp_A2000c30_n65_r12_v__",
        A1_sim          = [
            Ma_xp_A2000c30_n65_r12_v00,
            Ma_xp_A2000c30_n65_r12_v02,
            Ma_xp_A2000c30_n65_r12_v04,
            Ma_xp_A2000c30_n65_r12_v06,
            Ma_xp_A2000c30_n65_r12_v08,
            Ma_xp_A2000c30_n65_r12_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000c30_n65_r16_v__      = SimSet(
        name            = "set_Ma_xp_A2000c30_n65_r16_v__",
        A1_sim          = [
            Ma_xp_A2000c30_n65_r16_v00,
            Ma_xp_A2000c30_n65_r16_v02,
            Ma_xp_A2000c30_n65_r16_v04,
            Ma_xp_A2000c30_n65_r16_v06,
            Ma_xp_A2000c30_n65_r16_v08,
            Ma_xp_A2000c30_n65_r16_v10,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=v_inf"],
    )
    set_Ma_xp_A2000c30_n65_r12r16_v__   = SimSet(
        name            = "set_Ma_xp_A2000c30_n65_r12r16_v__",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000c30_n65_r12_v__.A1_sim,
            set_Ma_xp_A2000c30_n65_r16_v__.A1_sim,
        )),
    )

    # Resolutions
    set_Ma_xp_A2000_n___r12_v00     = SimSet(
        name            = "set_Ma_xp_A2000_n___r12_v00",
        A1_sim          = [
            Ma_xp_A2000_n50_r12_v00,
            Ma_xp_A2000_n55_r12_v00,
            Ma_xp_A2000_n60_r12_v00,
            Ma_xp_A2000_n65_r12_v00,
            Ma_xp_A2000_n70_r12_v00,
        ],
        A1_colour       = [cmap_rbow],
        A1_time         = [-1],
        A1_label        = [r"$10^{%.2g}$" % n for n in np.arange(5.0, 7.01, 0.5)],
    )
    set_Ma_xp_A2000_n___r14_v00     = SimSet(
        name            = "set_Ma_xp_A2000_n___r14_v00",
        copy            = "set_Ma_xp_A2000_n___r12_v00",
        A1_sim          = [
            Ma_xp_A2000_n50_r14_v00,
            Ma_xp_A2000_n55_r14_v00,
            Ma_xp_A2000_n60_r14_v00,
            Ma_xp_A2000_n65_r14_v00,
            Ma_xp_A2000_n70_r14_v00,
        ],
    )
    set_Ma_xp_A2000_n___r16_v00     = SimSet(
        name            = "set_Ma_xp_A2000_n___r16_v00",
        copy            = "set_Ma_xp_A2000_n___r12_v00",
        A1_sim          = [
            Ma_xp_A2000_n50_r16_v00,
            Ma_xp_A2000_n55_r16_v00,
            Ma_xp_A2000_n60_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2000_n70_r16_v00,
        ],
    )

    set_Ma_xp_A2000_n___r16_v00_r_  = SimSet(
        name            = "set_Ma_xp_A2000_n___r16_v00_r_",
        copy            = "set_Ma_xp_A2000_n___r16_v00",
        A1_sim          = [
            Ma_xp_A2000_n50_r16_v00,
            Ma_xp_A2000_n55_r16_v00,
            Ma_xp_A2000_n60_r16_v00,
            Ma_xp_A2000_n65_r16_v00,
            Ma_xp_A2000_n65_r16_v00_r1,
            Ma_xp_A2000_n65_r16_v00_r2,
            Ma_xp_A2000_n65_r16_v00_r3,
            Ma_xp_A2000_n65_r16_v00_r4,
            Ma_xp_A2000_n65_r16_v00_r5,
            Ma_xp_A2000_n65_r16_v00_r6,
            Ma_xp_A2000_n65_r16_v00_r7,
            Ma_xp_A2000_n70_r16_v00,
        ],
        A1_colour       = [
            cmap_rbow(np.linspace(0, 1, 5))[i] for i in [0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4]
        ],
        A1_label        = [
            r"$10^{%.2g}$" % n for n in [5, 5.5, 6, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 7]
        ],
    )

    set_Ma_xp_A2000_n___r___v00     = SimSet(
        name            = "set_Ma_xp_A2000_n___r___v00",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_n___r12_v00.A1_sim,
            set_Ma_xp_A2000_n___r14_v00.A1_sim,
            set_Ma_xp_A2000_n___r16_v00.A1_sim,
        )),
    )

    # Paper combinations etc
    set_Ma_xp_A2000_s___z_n65_r___v__   = SimSet(
        name            = "set_Ma_xp_A2000_s___z_n65_r___v__",
        A1_sim          = np.concatenate((
            set_Ma_xp_A2000_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_s047z_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_n65_r___v02.A1_sim,
            set_Ma_xp_A2000_n65_r___v04.A1_sim,
            set_Ma_xp_A2000_n65_r___v06.A1_sim,
            set_Ma_xp_A2000_n65_r___v08.A1_sim,
            set_Ma_xp_A2000c30_n65_r___v00.A1_sim,
        )),
        Di_param_Di_p_colour    = {
            "v_inf" : {
                np.round(0 + i * 0.2, 1): cmap_rbow(np.linspace(0, 1, 5))[i] for i in range(5)
            },
        },
        Di_param_Di_p_linestyle = {
            "L_2"   : {0: "-", 0.5: ls_dash, 1: ls_dash_dot, -1: ls_dash_dot_dot},
        },
    )

    # All paper results
    set_phodei_tables   = SimSet(
        name        = "set_phodei_tables",
        A1_sim      = np.concatenate((
            set_Ma_xp_A2000_n65_r___v__.A1_sim,
            set_Ma_xp_A2000_n___r___v00.A1_sim,
            set_Ma_xp_A_____n___r12r16_v00.A1_sim,
            set_Ma_xp_A2000c30_n65_r___v00.A1_sim,
            set_Ma_xp_A2000c30_n65_r12r16_v__.A1_sim,
            set_Ma_xp_A2000_s___x_n65_r12_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___y_n65_r12_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___mz_n65_r12_v00.A1_sim,
            set_Ma_xp_A2000_s___z_n65_r12_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___x_n65_r16_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___y_n65_r16_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___mz_n65_r16_v00.A1_sim,
            set_Ma_xp_A2000_s___z_n65_r16_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___x_n65_r20_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___y_n65_r20_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___z_n65_r20_v00.A1_sim[::-1],
            set_Ma_xp_A2000_s___x_n65_r16_v06.A1_sim[::-1],
            set_Ma_xp_A2000_s___y_n65_r16_v06.A1_sim[::-1],
            set_Ma_xp_A2000_s___z_n65_r16_v06.A1_sim[::-1],
            set_Ma_xp_A2000_s047z_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r___v00.A1_sim,
            set_Ma_xp_A2000_s047z_n65_r12_v__.A1_sim,
            set_Ma_xp_A2000_s047z_n65_r16_v__.A1_sim,
            set_Ma_xp_A2000_s047z_n65_r20_v__.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r12_v__.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r16_v__.A1_sim,
            set_Ma_xp_A2000_s030z_n65_r20_v__.A1_sim,
        )),
    )

# Rebound simulations
if True:
    # Base
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        po                  = plot_opt_phodei_reb,
        dir_proj            = ut.dir_gihr + "rebound/phodei/",
        category            = "phodei",
        t_end               = 5000,  # yr
        t_out_step          = 0.01,  # yr
        t_wall_restart      = 10,  # h
        t_2                 = 10,  # yr
        t_out_step_2        = 0.1,  # yr
        t_3                 = 500,  # yr
        t_out_step_3        = 0.5,  # yr
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, A1_pos=[Ma.o.a, 0, 0], A1_vel=[0, Ma.o.v_mean, 0], name="Mars"),
        ],
        id_cent             = 1,
        id_orb_start        = 2,
        sim_in              = Ma_xp_A2000_n65_r16_v00,
        time_in             = -1,
        # m_fof_min           = 1e16,  # kg (--> 76)
        m_fof_min           = 3e15,  # kg (--> 105)
        bnd_orb_Q_max       = 350 * Ma.R,  # 1.207 R_Hill
        test_masses         = "all",
        sim_in_rot_ax       = "x",
        sim_in_angle        = -0,   # inclination, i (negative)
        sim_in_rot_ax_2     = "z",
        sim_in_angle_2      = 0,   # orientation, phi
        oblate_id           = 1,
        oblate_J2           = 1960.45e-6,
        oblate_obliquity    = -30,  # Mars obliquity, o (negative)
        id_rel_prim         = 1,
        rel_distance_max    = 2 * Ma.o.R_Hill,
        collision_mode      = ur.reb_collision_mode_print,
        t_collision_delay   = 0.02,
        id_frame_shift      = 1,
    )

    # Examples and tests
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        oblate_obliquity    = -30,  # o (negative)
        # Speed tests: wall time to 10 yr (min), estimated time to 10 kyr (h)
        # t_end               = 10,  # yr
        # num_fof             = 30,   # 1.9, 32
        # num_fof             = 40,   # 3.0, 50
        # num_fof             = 50,   # 4.3, 72
        # num_fof             = 60,   # 5.6, 93
        # num_fof             = 70,   # 7.6, 126
        # num_fof             = 76,   # 8.5, 142
        # num_fof             = 80,   # 9.2, 153
        # num_fof             = 90,   # 11.4, 190
        # num_fof             = 100,  # 13.7, 228
        # num_fof             = 110,  # 16.2, 270
    )

    # Resolutions
    reb_Ma_xp_A2000_n70_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_n70_r12_v00,
        m_fof_min           = 3e15,  # kg (--> 191)
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r12_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r12_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r12_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 75,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r12_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )

    reb_Ma_xp_A2000_n70_r14_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r14_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_n70_r14_v00,
        m_fof_min           = 3e15,  # kg (--> 198)
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r14_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r14_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r14_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 75,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r14_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r14_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r14_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )

    reb_Ma_xp_A2000_n70_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_n70_r16_v00,
        m_fof_min           = 3e15,  # kg (--> 164)
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 75,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n70_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n70_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n70_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        oblate_obliquity    = -30,  # o (negative)
        Di_misc             = {
            "A1_idx_orb"    : [0, 1, 2, 4, 5, 6, 8, 10, 13, 15, 17, 18, 22, 23, 28],
        }
    )

    # Encounter periapsis angles, no inclination
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 15,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s030z
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s030z
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s030x
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s030x
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s030y
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s030y
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s047z
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s047z
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s047mz
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s047mz_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s047mz
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s047x
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s047x
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, no inclination, s047y
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r16_v00,
        sim_in_angle        = -0,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i00_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i00_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Encounter periapsis angles, mid inclination, s047y
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        oblate_obliquity    = -30,  # o (negative)
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 15,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 45,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 75,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 105,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 120,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 135,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 150,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 165,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 180,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 195,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 210,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 225,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 240,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 255,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 270,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 285,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 300,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 315,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 330,  # phi
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i15_o30",
        copy                = "reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30",
        sim_in_angle_2      = 345,  # phi
    )

    # Inclinations, for a few periapsis angles
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i05_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i05_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 0,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i30_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i30_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i45_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i45_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i60_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i60_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i75_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i75_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i90_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i90_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi030_i05_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i05_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i30_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i30_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i45_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i45_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i60_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i60_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i75_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i75_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i90_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i90_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi060_i05_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i05_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i30_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i30_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i45_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i45_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i60_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i60_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i75_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i75_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i90_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i90_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi090_i05_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i05_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i30_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i30_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i45_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i45_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i60_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i60_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i75_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i75_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i90_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i90_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    # Inclinations, for a few periapsis angles, s047z
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i05_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i05_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 0,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i30_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i30_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i45_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i45_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i60_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i60_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i75_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i75_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i90_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i90_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i05_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i05_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 30,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i30_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i30_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i45_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i45_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i60_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i60_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i75_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i75_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i90_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i90_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i05_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i05_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 60,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i30_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i30_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i45_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i45_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i60_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i60_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i75_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i75_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i90_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i90_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i05_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i05_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -5,   # i (negative)
        sim_in_angle_2      = 90,   # phi
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i30_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i30_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -30,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i45_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i45_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -45,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i60_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i60_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -60,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i75_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i75_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -75,   # i (negative)
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i90_o30    = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i90_o30",
        copy                = "reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30",
        sim_in_angle        = -90,   # i (negative)
    )

    # Disruption scenario: periapsis distance, for a few periapsis angles
    reb_Ma_xp_A2000_n65_r11_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r11_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A2000_n65_r11_v00,
    )
    reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A2000_n65_r12_v00,
    )
    reb_Ma_xp_A2000_n65_r13_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r13_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A2000_n65_r13_v00,
    )
    reb_Ma_xp_A2000_n65_r14_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r14_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r14_v00,
    )
    reb_Ma_xp_A2000_n65_r15_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r15_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A2000_n65_r15_v00,
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_n65_r17_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r17_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r17_v00,
    )
    reb_Ma_xp_A2000_n65_r18_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r18_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r18_v00,
    )
    reb_Ma_xp_A2000_n65_r19_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r19_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r19_v00,
    )

    reb_Ma_xp_A2000_n65_r11_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r11_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A2000_n65_r11_v00,
    )
    reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A2000_n65_r12_v00,
    )
    reb_Ma_xp_A2000_n65_r13_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r13_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A2000_n65_r13_v00,
    )
    reb_Ma_xp_A2000_n65_r14_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r14_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r14_v00,
    )
    reb_Ma_xp_A2000_n65_r15_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r15_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A2000_n65_r15_v00,
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_n65_r17_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r17_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r17_v00,
    )
    reb_Ma_xp_A2000_n65_r18_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r18_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r18_v00,
    )
    reb_Ma_xp_A2000_n65_r19_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r19_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r19_v00,
    )

    reb_Ma_xp_A2000_n65_r11_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r11_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A2000_n65_r11_v00,
    )
    reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A2000_n65_r12_v00,
    )
    reb_Ma_xp_A2000_n65_r13_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r13_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A2000_n65_r13_v00,
    )
    reb_Ma_xp_A2000_n65_r14_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r14_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r14_v00,
    )
    reb_Ma_xp_A2000_n65_r15_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r15_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A2000_n65_r15_v00,
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r17_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r17_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r17_v00,
    )
    reb_Ma_xp_A2000_n65_r18_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r18_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r18_v00,
    )
    reb_Ma_xp_A2000_n65_r19_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r19_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r19_v00,
    )

    reb_Ma_xp_A2000_n65_r11_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r11_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        sim_in              = Ma_xp_A2000_n65_r11_v00,
    )
    reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        sim_in              = Ma_xp_A2000_n65_r12_v00,
    )
    reb_Ma_xp_A2000_n65_r13_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r13_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        sim_in              = Ma_xp_A2000_n65_r13_v00,
    )
    reb_Ma_xp_A2000_n65_r14_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r14_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r14_v00,
    )
    reb_Ma_xp_A2000_n65_r15_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r15_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r15_v00,
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r17_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r17_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r17_v00,
    )
    reb_Ma_xp_A2000_n65_r18_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r18_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r18_v00,
    )
    reb_Ma_xp_A2000_n65_r19_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r19_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r19_v00,
    )

    # Disruption scenario: speed, for a few periapsis angles
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_n65_r16_v02_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v02_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v02,
    )
    reb_Ma_xp_A2000_n65_r16_v04_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v04_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v04,
    )
    reb_Ma_xp_A2000_n65_r16_v06_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v06_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v06,
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_n65_r16_v02_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v02_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v02,
    )
    reb_Ma_xp_A2000_n65_r16_v04_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v04_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v04,
    )
    reb_Ma_xp_A2000_n65_r16_v06_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v06_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v06,
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r16_v02_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v02_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v02,
    )
    reb_Ma_xp_A2000_n65_r16_v04_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v04_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v04,
    )
    reb_Ma_xp_A2000_n65_r16_v06_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v06_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v06,
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r16_v02_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v02_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v02,
    )
    reb_Ma_xp_A2000_n65_r16_v04_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v04_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v04,
    )
    reb_Ma_xp_A2000_n65_r16_v06_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v06_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v06,
    )

    # Disruption scenario: spin, for a few periapsis angles
    #   r12
    #       phi000
    reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r12_v00,
    )

    #       phi030
    reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r12_v00,
    )

    #       phi060
    reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r12_v00,
    )

    #       phi090
    reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r12_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r12_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r12_v00,
    )

    #   r16
    #       phi000
    reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s036z_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s086z_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi000_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s036x_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s086x_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s036y_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s086y_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r16_v00,
    )

    #       phi030
    reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s036z_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s086z_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi030_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s036x_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s086x_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s036y_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s086y_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r16_v00,
    )

    #       phi060
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s036z_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s086z_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi060_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s036x_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s086x_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s036y_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s086y_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r16_v00,
    )

    #       phi090
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s036z_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s086z_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086mz_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi090_i15_o30 = RebSim(
        name                = "reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170mz_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s036x_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s086x_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r16_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s036y_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s086y_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r16_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r16_v00,
    )

    #   r20 (non-spin objects just for copying, not to run since no initial capture)
    reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r20_v00,
    )
    reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r20_v00,
    )
    reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r20_v00,
    )
    reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r20_v00,
    )
    #       phi000
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r20_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r20_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r20_v00,
    )

    #       phi030
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r20_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r20_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r20_v00,
    )

    #       phi060
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r20_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r20_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r20_v00,
    )

    #       phi090
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036z_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036z_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086z_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086z_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086z_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170z_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170z_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170z_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030x_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030x_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036x_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036x_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047x_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047x_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086x_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086x_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086x_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170x_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170x_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170x_n65_r20_v00,
    )

    reb_Ma_xp_A2000_s030y_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030y_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s036y_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s036y_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s036y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s047y_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047y_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s086y_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s086y_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s086y_n65_r20_v00,
    )
    reb_Ma_xp_A2000_s170y_n65_r20_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s170y_n65_r20_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s170y_n65_r20_v00,
    )

    # Disruption scenario: spin periapsis distance, for a few periapsis angles
    #   s030z
    reb_Ma_xp_A2000_s030z_n65_r11_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r11_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r13_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r13_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r14_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r14_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r15_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r15_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r17_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r17_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r18_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r18_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r19_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r19_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r22_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r22_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r24_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r24_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s030z_n65_r11_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r11_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r13_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r13_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r14_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r14_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r15_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r15_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r17_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r17_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r18_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r18_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r19_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r19_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r22_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r22_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r24_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r24_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s030z_n65_r11_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r11_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r13_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r13_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r14_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r14_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r15_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r15_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r17_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r17_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r18_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r18_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r19_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r19_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r22_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r22_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r24_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r24_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s030z_n65_r11_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r11_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r12_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r13_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r13_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r14_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r14_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r15_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r15_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r17_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r17_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r18_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r18_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r19_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r19_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r20_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s030z_n65_r22_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r22_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s030z_n65_r24_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s030z_n65_r24_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s030z_n65_r24_v00,
    )

    #   s047z
    reb_Ma_xp_A2000_s047z_n65_r11_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r11_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r13_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r13_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r14_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r14_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r15_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r15_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r17_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r17_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r18_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r18_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r19_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r19_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi000_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r22_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r22_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r24_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r24_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s047z_n65_r11_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r11_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r13_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r13_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r14_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r14_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r15_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r15_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r17_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r17_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r18_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r18_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r19_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r19_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi030_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r22_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r22_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r24_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r24_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s047z_n65_r11_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r11_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r13_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r13_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r14_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r14_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r15_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r15_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r17_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r17_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r18_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r18_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r19_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r19_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi060_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r22_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r22_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r24_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r24_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r24_v00,
    )

    reb_Ma_xp_A2000_s047z_n65_r11_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r11_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r11_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r12_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r13_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r13_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r13_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r14_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r14_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r14_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r15_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r15_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r15_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r17_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r17_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r17_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r18_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r18_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r18_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r19_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r19_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r19_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r20_v00_phi090_i15_o30
    reb_Ma_xp_A2000_s047z_n65_r22_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r22_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r22_v00,
    )
    reb_Ma_xp_A2000_s047z_n65_r24_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_s047z_n65_r24_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_s047z_n65_r24_v00,
    )

    # Mars obliquity
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        oblate_obliquity    = -0,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o15  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o15",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00",
        oblate_obliquity    = -15,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o45  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o45",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00",
        oblate_obliquity    = -45,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o60  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o60",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00",
        oblate_obliquity    = -60,  # o (negative)
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        oblate_obliquity    = -0,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o15  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o15",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00",
        oblate_obliquity    = -15,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o45  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o45",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00",
        oblate_obliquity    = -45,  # o (negative)
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o60  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o60",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00",
        oblate_obliquity    = -60,  # o (negative)
    )

    # Mars eccentricity
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu000   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu000",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=0 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu045   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu045",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=45 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu090   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu090",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=90 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu135   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu135",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=135 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu180   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu180",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=180 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu225   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu225",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=225 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu270   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu270",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=270 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu315   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu315",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=315 * deg_to_rad, name="Mars"),
        ],
    )

    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu000   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu000",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=0 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu045   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu045",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=45 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu090   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu090",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=90 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu135   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu135",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=135 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu180   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu180",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=180 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu225   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu225",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=225 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu270   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu270",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=270 * deg_to_rad, name="Mars"),
        ],
    )
    reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu315   = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu315",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        A1_picle            = [
            RebPicle(m=Su.M, R=Su.R, A1_pos=np.zeros(3), A1_vel=np.zeros(3), name="Sun"),
            RebPicle(m=Ma.M, R=Ma.R, a=1.524 * au_to_m, e=0.1, nu=315 * deg_to_rad, name="Mars"),
        ],
    )

    # Reoriented repeats
    reb_Ma_xp_A2000_n65_r16_v00_r1_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r1_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r1,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r2_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r2_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r2,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r3_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r3_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r3,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r4_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r4_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r4,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r5_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r5_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r5,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r6_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r6_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r6,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r7_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r7_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r7,
    )

    reb_Ma_xp_A2000_n65_r16_v00_r1_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r1_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r1,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r2_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r2_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r2,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r3_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r3_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r3,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r4_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r4_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r4,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r5_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r5_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r5,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r6_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r6_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r6,
    )
    reb_Ma_xp_A2000_n65_r16_v00_r7_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2000_n65_r16_v00_r7_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2000_n65_r16_v00_r7,
    )

    # Masses (fixed m_fof_min = 3e15)
    reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A1800_n60_r12_v00,  # (--> 14)
    )
    reb_Ma_xp_A1850_n60_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r12_v00,  # (--> 22)
    )
    reb_Ma_xp_A1900_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r12_v00,  # (--> 48)
    )
    reb_Ma_xp_A1950_n65_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r12_v00,  # (--> 86)
    )
    reb_Ma_xp_A2050_n70_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r12_v00,  # (--> 400)
    )
    reb_Ma_xp_A2100_n70_r12_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r12_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r12_v00,  # (--> 484)
    )

    reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A1800_n60_r12_v00,
    )
    reb_Ma_xp_A1850_n60_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r12_v00,
    )
    reb_Ma_xp_A1900_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r12_v00,
    )
    reb_Ma_xp_A1950_n65_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r12_v00,
    )
    reb_Ma_xp_A2050_n70_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r12_v00,
    )
    reb_Ma_xp_A2100_n70_r12_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r12_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r12_v00,
    )

    reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A1800_n60_r12_v00,
    )
    reb_Ma_xp_A1850_n60_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r12_v00,
    )
    reb_Ma_xp_A1900_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r12_v00,
    )
    reb_Ma_xp_A1950_n65_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r12_v00,
    )
    reb_Ma_xp_A2050_n70_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r12_v00,
    )
    reb_Ma_xp_A2100_n70_r12_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r12_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r12_v00,
    )

    reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        sim_in              = Ma_xp_A1800_n60_r12_v00,
    )
    reb_Ma_xp_A1850_n60_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r12_v00,
    )
    reb_Ma_xp_A1900_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r12_v00,
    )
    reb_Ma_xp_A1950_n65_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r12_v00,
    )
    reb_Ma_xp_A2050_n70_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r12_v00,
    )
    reb_Ma_xp_A2100_n70_r12_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r12_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r12_v00,
    )

    reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 0,   # phi
        sim_in              = Ma_xp_A1800_n60_r16_v00,  # (--> 3)
    )
    reb_Ma_xp_A1850_n60_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r16_v00,  # (--> 11)
    )
    reb_Ma_xp_A1900_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r16_v00,  # (--> 29)
    )
    reb_Ma_xp_A1950_n65_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r16_v00,  # (--> 56)
    )
    reb_Ma_xp_A2050_n70_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r16_v00,  # (--> 313)
    )
    reb_Ma_xp_A2100_n70_r16_v00_phi000_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r16_v00_phi000_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r16_v00,  # (--> 364)
    )

    reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 30,   # phi
        sim_in              = Ma_xp_A1800_n60_r16_v00,
    )
    reb_Ma_xp_A1850_n60_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r16_v00,
    )
    reb_Ma_xp_A1900_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r16_v00,
    )
    reb_Ma_xp_A1950_n65_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r16_v00,
    )
    reb_Ma_xp_A2050_n70_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r16_v00,
    )
    reb_Ma_xp_A2100_n70_r16_v00_phi030_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r16_v00_phi030_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r16_v00,
    )

    reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 60,   # phi
        sim_in              = Ma_xp_A1800_n60_r16_v00,
    )
    reb_Ma_xp_A1850_n60_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r16_v00,
    )
    reb_Ma_xp_A1900_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r16_v00,
    )
    reb_Ma_xp_A1950_n65_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r16_v00,
    )
    reb_Ma_xp_A2050_n70_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r16_v00,
    )
    reb_Ma_xp_A2100_n70_r16_v00_phi060_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r16_v00_phi060_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r16_v00,
    )

    reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30",
        sim_in_angle        = -15,   # i (negative)
        sim_in_angle_2      = 90,   # phi
        sim_in              = Ma_xp_A1800_n60_r16_v00,
    )
    reb_Ma_xp_A1850_n60_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1850_n60_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1850_n60_r16_v00,
    )
    reb_Ma_xp_A1900_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1900_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1900_n65_r16_v00,
    )
    reb_Ma_xp_A1950_n65_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A1950_n65_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A1950_n65_r16_v00,
    )
    reb_Ma_xp_A2050_n70_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2050_n70_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2050_n70_r16_v00,
    )
    reb_Ma_xp_A2100_n70_r16_v00_phi090_i15_o30  = RebSim(
        name                = "reb_Ma_xp_A2100_n70_r16_v00_phi090_i15_o30",
        copy                = "reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30",
        sim_in              = Ma_xp_A2100_n70_r16_v00,
    )

# Rebound simulation sets
if True:
    # Encounter periapsis angles, without and with spin
    set_reb_Ma_xp_A2000_n65_r16_v00_phi____i00_o30      = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi____i15_o30      = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi____i00_o30_2    = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi____i00_o30_2",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi315_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=phi"],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi____i15_o30_2    = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi____i15_o30_2",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi315_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = ["auto=phi"],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s030z_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s030z_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s030x_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030x_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s030x_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030x_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s030y_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030y_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s030y_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030y_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi____i00_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi____i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s047x_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047x_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047x_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047x_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s047y_n65_r16_v00_phi____i00_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047y_n65_r16_v00_phi____i00_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i00_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i00_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "0",
        set_label_title = r"$i$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047y_n65_r16_v00_phi____i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047y_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi105_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi120_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi135_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi150_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi165_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi180_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi195_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi210_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi225_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi240_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi255_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi270_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi285_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi300_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi315_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi330_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi345_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        set_label       = "15",
        set_label_title = r"$i$ ($^\circ$)",
    )

    # Inclinations, for a few periapsis angles
    set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i05_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i45_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i75_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i90_o30,
        ],
        A1_label        = [
            r"$0$", r"$5$", r"$15$", r"$30$", r"$45$", r"$60$", r"$75$", r"$90$",
            r"$105$", r"$120$", r"$135$", r"$150$", r"$165$", r"$180$",
        ],
        legend_title    = r"$i$ ($^\circ$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi030_i___o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi030_i___o30",
        copy            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i05_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i45_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i75_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i90_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i___o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i___o30",
        copy            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i05_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i45_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i75_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i90_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i___o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i___o30",
        copy            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i05_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i45_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i75_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i90_o30,
        ],
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i05_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i45_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i75_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i90_o30,
        ],
        A1_label        = [
            r"$0$", r"$5$", r"$15$", r"$30$", r"$45$", r"$60$", r"$75$", r"$90$",
            r"$105$", r"$120$", r"$135$", r"$150$", r"$165$", r"$180$",
        ],
        legend_title    = r"$i$ ($^\circ$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i___o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i05_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i45_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i75_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i90_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i___o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i05_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i45_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i75_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i90_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i___o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i05_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i45_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i75_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i90_o30,
        ],
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_n65_r16_v00_phi____i___o30          = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi____i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i90_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i90_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i90_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i30_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i60_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i90_o30,
        ],
        A1_colour       = np.repeat(
            [cmap_rbow(np.linspace(0, 1, 4)[i]) for i in range(4)], 4, axis=0
        ),
        A1_linestyle    = np.tile([ls_dot, ls_dash, ls_dash_dot, "-"], 4),
    )
    set_reb_Ma_xp_A2000__s047z_n65_r16_v00_phi____i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000__s047z_n65_r16_v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi015_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi045_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,
        ],
        A1_colour       = [cmap_rbow(np.linspace(0, 1, 7)[i]) for i in range(7)] * 2,
        A1_linestyle    = ["-"] * 7 + [ls_dash] * 7,
    )
    set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i___o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i90_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i90_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i90_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i60_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i90_o30,
        ],
        A1_colour       = np.repeat(
            [cmap_rbow(np.linspace(0, 1, 4)[i]) for i in range(4)], 4, axis=0
        ),
        A1_linestyle    = np.tile([ls_dot, ls_dash, ls_dash_dot, "-"], 4),
    )
    set_Ma_xp_A2000_s047__n65_r16_v00_                      = SimSet(
        name            = "set_Ma_xp_A2000_s047__n65_r16_v00_",
        A1_sim          = [
            Ma_xp_A2000_n70_r16_v00,
            Ma_xp_A2000_s047x_n65_r16_v00,
            Ma_xp_A2000_s047y_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,

            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
            Ma_xp_A2000_s047z_n65_r16_v00,
            # Ma_xp_A2000_s047z_n65_r16_v00,
            # Ma_xp_A2000_s047z_n65_r16_v00,
        ],
        A1_colour       = ["k", A1_c[0], A1_c[1], A1_c[2]] + [A1_c[2]] * 4,
        A1_label        = [r"$0$", r"$x$", r"$y$", r"$z$"] + [None] * 4,
        A1_linestyle    = ["-"] * 4 + [
            ls_dash,
            (0, (7, 2, 1, 2)),
            (0, (7, 2, 1, 2, 1, 2)),
            (0, (7, 2, 1, 2, 1, 2, 1, 2)),
            # (0, (7, 2, 1, 2, 1, 2, 1, 2, 1, 2)),
            # (0, (7, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2)),
        ],
        legend_title    = r"Spin Axis",
    )
    set_reb_Ma_xp_A2000_s047__n65_r16_v00_phi____i___o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047__n65_r16_v00_phi____i___o30",
        copy            = "set_Ma_xp_A2000_s047__n65_r16_v00_",
        A1_sim          = [
            reb_Ma_xp_A2000_n70_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,

            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i00_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i30_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i45_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i60_o30,
            # reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i75_o30,
            # reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i90_o30,
        ],
    )

    # Disruption scenario: periapsis distance and speed
    set_reb_Ma_xp_A2000_n65_r___v00_phi000_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r11_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r13_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r14_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r15_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r17_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r18_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r19_v00_phi000_i15_o30,
        ],
        A1_label        = [r"$%.1f$" % r for r in np.arange(1.1, 1.91, 0.1)],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r___v00_phi030_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r___v00_phi030_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r11_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r13_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r14_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r15_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r17_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r18_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r19_v00_phi030_i15_o30,
        ],
        A1_label        = [r"$%.1f$" % r for r in np.arange(1.1, 1.91, 0.1)],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_n65_r___v00_phi060_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r___v00_phi060_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r11_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r13_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r14_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r15_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r17_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r18_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r19_v00_phi060_i15_o30,
        ],
        A1_label        = [r"$%.1f$" % r for r in np.arange(1.1, 1.91, 0.1)],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_n65_r___v00_phi090_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r___v00_phi090_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r11_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r13_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r14_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r15_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r17_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r18_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r19_v00_phi090_i15_o30,
        ],
        A1_label        = [r"$%.1f$" % r for r in np.arange(1.1, 1.91, 0.1)],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    set_reb_Ma_xp_A2000_n65_r16_v___phi000_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v___phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v02_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v04_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v06_phi000_i15_o30,
        ],
        A1_label        = [r"$0$", r"$0.2$", r"$0.4$", r"$0.6$"],
        legend_title    = r"$v_{\infty}$ (km s$^{-1}$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v___phi030_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v___phi030_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v02_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v04_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v06_phi030_i15_o30,
        ],
        A1_label        = [r"$0$", r"$0.2$", r"$0.4$", r"$0.6$"],
        legend_title    = r"$v_{\infty}$ (km s$^{-1}$)",
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_n65_r16_v___phi060_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v___phi060_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v02_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v04_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v06_phi060_i15_o30,
        ],
        A1_label        = [r"$0$", r"$0.2$", r"$0.4$", r"$0.6$"],
        legend_title    = r"$v_{\infty}$ (km s$^{-1}$)",
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_n65_r16_v___phi090_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v___phi090_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v02_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v04_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v06_phi090_i15_o30,
        ],
        A1_label        = [r"$0$", r"$0.2$", r"$0.4$", r"$0.6$"],
        legend_title    = r"$v_{\infty}$ (km s$^{-1}$)",
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    # Disruption scenario: spin
    set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi000_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi030_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi060_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi090_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030mz_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036mz_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047mz_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi000_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi030_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi060_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi090_i15_o30   = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086mz_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170mz_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{z}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___x_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030x_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036x_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047x_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086x_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170x_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi000_i15_o30,
        ],
        A1_label        = ["1", "3/4", "1/2", "1/4", "1/8", "0"],  # L_max
        legend_title    = r"$L_{x}$ $(L_{\rm max})$",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s___y_n65_r20_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030y_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s036y_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047y_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s086y_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s170y_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r20_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    # Disruption scenario: z-spin periapsis distance
    set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r11_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r13_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r14_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r15_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r17_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r18_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r19_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r22_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r24_v00_phi000_i15_o30,
        ],
        A1_label        = [
            r"$%.1f$" % r for r in np.append(np.arange(1.1, 2.01, 0.1), [2.2, 2.4])
        ],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r11_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r13_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r14_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r15_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r17_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r18_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r19_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r22_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r24_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r11_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r13_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r14_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r15_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r17_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r18_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r19_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r22_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r24_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s030z_n65_r11_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r13_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r14_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r15_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r17_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r18_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r19_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r22_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s030z_n65_r24_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r11_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r13_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r14_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r15_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r17_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r18_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r19_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r22_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r24_v00_phi000_i15_o30,
        ],
        A1_label        = [
            r"$%.1f$" % r for r in np.append(np.arange(1.1, 2.01, 0.1), [2.2, 2.4])
        ],
        legend_title    = r"$q$ ($R_{\mars{}\!}$)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi030_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r11_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r13_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r14_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r15_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r17_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r18_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r19_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r22_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r24_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi060_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r11_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r13_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r14_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r15_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r17_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r18_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r19_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r22_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r24_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi090_i15_o30    = SimSet(
        name            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_s047z_n65_r11_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r13_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r14_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r15_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r17_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r18_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r19_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r20_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r22_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_s047z_n65_r24_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    # Mars obliquity
    set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o__  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o__",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o00,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o15,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o45,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o60,
        ],
        A1_label        = [r"$0$", r"$15$", r"$30$", r"$45$", r"$60$"],
        legend_title    = r"$o_{\mars{}\!}$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o__  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o__",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o00,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o15,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o45,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o60,
        ],
        A1_label        = [r"$0$", r"$15$", r"$30$", r"$45$", r"$60$"],
        legend_title    = r"$o_{\mars{}\!}$ ($^\circ$)",
    )

    # Mars eccentricity
    set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu___  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu___",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu000,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu045,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu090,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu135,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu180,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu225,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu270,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30_e010_nu315,
        ],
        A1_label        = [r"N/A"] + [r"$%d$" % nu for nu in np.arange(0, 315.1, 45)],
        legend_title    = r"$\nu_{\mars{}\!}$ ($^\circ$)",
        set_label       = "60",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu___  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu___",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu000,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu045,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu090,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu135,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu180,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu225,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu270,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu315,
        ],
        A1_label        = [r"N/A"] + [r"$%d$" % nu for nu in np.arange(0, 315.1, 45)],
        legend_title    = r"$\nu_{\mars{}\!}$ ($^\circ$)",
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    # High resolutions
    set_reb_Ma_xp_A2000_n70_r___v00_phi____i15_o30      = SimSet(
        name            = "set_reb_Ma_xp_A2000_n70_r___v00_phi____i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n70_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n70_r12_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_n70_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n70_r14_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n70_r14_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_n70_r14_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n70_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n70_r16_v00_phi075_i15_o30,
            reb_Ma_xp_A2000_n70_r16_v00_phi090_i15_o30,
        ],
        A1_colour       = [
            cmap_rbow(np.linspace(0, 1, 3))[i] for i in [0, 0, 0, 1, 1, 1, 2, 2, 2]
        ],
        A1_linestyle    = [ls_dot, ls_dash, "-"] * 3,
    )

    # Reoriented repeats
    set_reb_Ma_xp_A2000_n65_r16_v00_r__phi060_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_r__phi060_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r1_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r2_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r3_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r4_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r5_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r6_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r7_phi060_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = [r"%d" % i for i in np.arange(8)],
        set_label       = "60",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A2000_n65_r16_v00_r__phi090_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A2000_n65_r16_v00_r__phi090_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r1_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r2_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r3_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r4_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r5_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r6_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_r7_phi090_i15_o30,
        ],
        A1_colour       = [cmap_rbow],
        A1_label        = [r"%d" % i for i in np.arange(8)],
        set_label       = "90",
        set_label_title = r"$\phi$ ($^\circ$)",
    )

    # Masses
    set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A1850_n60_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A1900_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A1950_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2050_n70_r12_v00_phi000_i15_o30,
            reb_Ma_xp_A2100_n70_r12_v00_phi000_i15_o30,
        ],
        A1_label        = [r"$10^{%g}$" for logm in [18, 18.5, 19, 19.5, 20, 20.5, 21]],
        legend_title    = r"$M_{\rm A}$ (kg)",
        set_label       = "0",
        set_label_title = r"$\phi$ ($^\circ$)",
    )
    set_reb_Ma_xp_A_____n___r12_v00_phi030_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r12_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A1850_n60_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A1900_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A1950_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2050_n70_r12_v00_phi030_i15_o30,
            reb_Ma_xp_A2100_n70_r12_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A_____n___r12_v00_phi060_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r12_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A1850_n60_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A1900_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A1950_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2050_n70_r12_v00_phi060_i15_o30,
            reb_Ma_xp_A2100_n70_r12_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A_____n___r12_v00_phi090_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r12_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A1850_n60_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A1900_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A1950_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2050_n70_r12_v00_phi090_i15_o30,
            reb_Ma_xp_A2100_n70_r12_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )
    set_reb_Ma_xp_A_____n___r16_v00_phi000_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r16_v00_phi000_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A1850_n60_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A1900_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A1950_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2050_n70_r16_v00_phi000_i15_o30,
            reb_Ma_xp_A2100_n70_r16_v00_phi000_i15_o30,
        ],
        set_label       = "0",
    )
    set_reb_Ma_xp_A_____n___r16_v00_phi030_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r16_v00_phi030_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A1850_n60_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A1900_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A1950_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2050_n70_r16_v00_phi030_i15_o30,
            reb_Ma_xp_A2100_n70_r16_v00_phi030_i15_o30,
        ],
        set_label       = "30",
    )
    set_reb_Ma_xp_A_____n___r16_v00_phi060_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r16_v00_phi060_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A1850_n60_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A1900_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A1950_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2050_n70_r16_v00_phi060_i15_o30,
            reb_Ma_xp_A2100_n70_r16_v00_phi060_i15_o30,
        ],
        set_label       = "60",
    )
    set_reb_Ma_xp_A_____n___r16_v00_phi090_i15_o30  = SimSet(
        name            = "set_reb_Ma_xp_A_____n___r16_v00_phi090_i15_o30",
        copy            = "set_reb_Ma_xp_A_____n___r12_v00_phi000_i15_o30",
        A1_sim          = [
            reb_Ma_xp_A1800_n60_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A1850_n60_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A1900_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A1950_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2050_n70_r16_v00_phi090_i15_o30,
            reb_Ma_xp_A2100_n70_r16_v00_phi090_i15_o30,
        ],
        set_label       = "90",
    )

    # All paper results
    set_phodei_reb_tables   = SimSet(
        name        = "set_phodei_reb_tables",
        A1_sim      = np.concatenate((
            set_reb_Ma_xp_A2000_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi000_i___o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi030_i___o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi060_i___o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i___o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r___v00_phi000_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r___v00_phi030_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r___v00_phi060_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r___v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v___phi000_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v___phi030_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v___phi060_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v___phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o__.A1_sim,
            set_reb_Ma_xp_A2000_n65_r16_v00_phi090_i15_o30_e010_nu___.A1_sim,
            set_reb_Ma_xp_A_____n___r12_v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A_____n___r16_v00_phi000_i15_o30.A1_sim,
            set_reb_Ma_xp_A_____n___r16_v00_phi030_i15_o30.A1_sim,
            set_reb_Ma_xp_A_____n___r16_v00_phi060_i15_o30.A1_sim,
            set_reb_Ma_xp_A_____n___r16_v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030x_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030y_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030z_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047x_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047y_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047mz_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047z_n65_r16_v00_phi____i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi000_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi030_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi060_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s030z_n65_r___v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi000_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi030_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi060_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s047z_n65_r___v00_phi090_i15_o30.A1_sim,
            # set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___x_n65_r12_v00_phi090_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___y_n65_r12_v00_phi090_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___mz_n65_r12_v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s___z_n65_r12_v00_phi090_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___x_n65_r16_v00_phi090_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___y_n65_r16_v00_phi090_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi000_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi030_i15_o30.A1_sim[::-1],
            # set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi060_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___mz_n65_r16_v00_phi090_i15_o30.A1_sim,
            set_reb_Ma_xp_A2000_s___z_n65_r16_v00_phi090_i15_o30.A1_sim[::-1],
            set_reb_Ma_xp_A2000_s___z_n65_r20_v00_phi090_i15_o30.A1_sim[:-1][::-1],
        )),
    )

# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR objects_phodei.py  ====\n")
