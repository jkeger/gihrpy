"""
Classes for organising GIHR projects, etc.

Objects of these classes are effectively used in place of input files.

See also classes_sim.py.

Copy
----
Most classes have an optional copy argument. This can be used to copy all the
attributes from an existing instance to this new one, if they are not set
directly for the new object. This makes it easy to create many similar versions
of objects with only small changes between each one.

See copy_object() in jkeger.py for more details.

Dictionaries
------------
Most classes automatically fill a dictionary of their instances each time a new
class object is made, with the object's name attribute as the key. This allows
functions that take one of these class objects as an argument to either accept
the object itself or the object's name to extract it from the dictionary.

Init
----
Most classes use __init__() only to set the input attributes, and have a
separate init() method that finishes setting up the object before use, to avoid
wasting time preparing objects that aren't needed.

Notation etc
------------
+ See README.md.
+ param_* is usually a parameter name (str), while p_* is a value (e.g. float).
"""

from jkeger import *
import gihrpy.utilities as ut

# Reference defaults
Di_param_unit_label_def = Di_param_unit_label


class PlotOptions:
    """Class for analysis and plotting options.

    Parameters
    ----------
    copy : PlotOptions
        A different PlotOptions object from which to copy all
        attributes that are not provided as other arguments here.

    category : str
        The base project type.

    dir_save : str
        The base folder for saving plots.

    Di_mat_colour : {str : str}
        A dictionary of material names and colours. See materials.py for
        the details and defaults. Is converted to Di_id_colour for
        usually easier use.

    Di_mat_highlt : {str : str}
        Like Di_mat_colour but for highlighting colours. Is converted
        to Di_id_highlight for usually easier use. Defaults to all
        c_highlight.

    Di_param_unit_label : {str : ParamUnitLabel}
        A dictionary of parameter unit/label/etc objects.

    unit_label_r_def, unit_label_m_def : ParamUnitLabel
        Unit/label/etc defaults for mass and distance parameters.

    Di_param_log : {str : bool}
        A dictionary of particle parameters and whether or not to plot with a
        logarithmic scale.

    Di_param_lim : {str : [float, float]}
        A dictionary of particle parameters and minimum, maximum values
        for setting axis limits on plots (arbitrary plot units).

    Di_misc_colour : {str : str}
        A dictionary of miscellanous colours with custom keys.

    dpi : int
        The dots-per-inch resolution for saving figures.

    fontsize : int
        The base fontsize, e.g. for axis labels.

    fontsize_text : int
        The fontsize for text annotations.

    cmap : mpl.Colormap
        A default matplotlib colour map for plotting.

    func_cmap_param : function
        A function with arguments of (po, param, A1_param) that returns
        (cmap, vmin, vmax, norm, extend).

    Di_cmap_param : {}
        A dictionary with any parameters required for func_cmap_param.

    A1_ax_lim : [float]
        The axis limits for a plot: [x_0, x_1, y_0, y_1] (r_unit), or None.

    tick_base, tick_base_minor : float
        The base separation for major and minor axis ticks (r_unit).

    func_ax_lim : function
        A function with arguments of (po, time, A2_pos) that returns custom-set
        (A1_ax_lim, tick_base, tick_base_minor).

    Di_ax_lim : {}
        A dictionary with any parameters required for func_ax_lim.

    thick : float
        The thickness of a slice in z about z=0 for plotting particles (r_unit).
        None    Plot all particles.
        -ve     Plot all particles between z_min and z=0 in |thick| slices.
        0       Plot the particles within their smoothing lengths of z=0.

    marker_size : float
        The size of the points to plot.

    func_marker_size : function
        A function with arguments of (po, A1_ax_lim) that returns a custom
        marker size.

    Di_marker_size : {}
        A dictionary with any parameters required for func_marker_size.

    ref_m_picle : float
        A reference particle mass. If the particle mass is different, e.g. for
        simulations at other resolutions, then scale things like the marker
        size accordingly.

    alpha : float
        The opacity for plotting particles.

    proj : str
        The plotting projection plane: xy, xz, zy, or yz.

    rgba_bg : [float]
        The background colour as rgba.

    text_colour : str
        The text colour.

    text_ul, text_uc, text_ur, text_ll, text_lc, text_lr : str
        Text for the upper/lower left/centre/right corners. Automatic options:
        tstamp          Time stamp (h).
        eval=*          An evaluatable string, e.g. latex format of a property.
        scale=*         Scale bar, length=* (r_unit). Set *=auto for automatic.

    Di_init_orbit : {}
        A dictionary with parameters required for plot_init_orbit.

    A1_po_inset : [PlotOptions]
        A list of other PlotOptions objects for an inset panel. Note: this is
        perhaps best set after initialising the parent object, so that the inset
        objects can set the parent as their copy.

    A1_inset_loc : [float]
        The location if this is an inset panel: [left, bottom, width, height],
        as fractions of the outer axes' width and height. Set width or height to
        None for an automatic aspect ratio to match the axis limits.

    fof_id_max : int
        The maximum friends of friends group ID to plot individually.
        Smaller groups will be accumulated together.

    fof_m_min : float
        The minimum friends of friends group mass (kg) to plot individually.
        Smaller groups will be accumulated together.

    fof_n_min : int
        The minimum number of particles required to be put in a group.

    R_Roche : float
        The Roche limit radius (m).

    R_Hill : float
        The Hill sphere radius (m).

    xz_size, zy_size : str
        The fractional size of a z-projection panels, e.g. "50%".

    thick_xz : float
        As above, for the z-projection panel.

    Fp_suffix : str
        An extra string to add to the filename when saving a plot, e.g. "_test".

    func_snapshot_add_misc : function
        A function with arguments of (po, time, param_c, ax) to add any
        miscellaneous elements to a snapshot plot.

    func_particle_params_extra : function
        A function with arguments of (sim, time, param_x, param_y, ax) to add
        any miscellaneous elements to a particle parameters plot.

    func_param_subsets_extra : function
        A function with arguments of (sim_set, param_x, param_y, param_c,
        param_ls, ax) to add any miscellaneous elements to a parameter subsets
        plot.

    func_reb_orbit_evol_extra : function
        A function with arguments of (rsim, A1_param_y, time_split,
        time_split_2, A1_ax_t0, A1_ax_t1, A1_ax_t2) to add any miscellaneous
        elements to a rebound orbit evolution plot.
    """

    def __init__(
        self,
        copy=None,
        category=None,
        dir_save=None,
        Di_mat_colour=None,
        Di_mat_highlt=None,
        Di_param_unit_label=None,
        unit_label_r_def=None,
        unit_label_m_def=None,
        Di_param_log=None,
        Di_param_lim=None,
        Di_misc_colour=None,
        dpi=None,
        fontsize=None,
        fontsize_text=None,
        cmap=None,
        func_cmap_param=None,
        Di_cmap_param=None,
        A1_ax_lim=None,
        tick_base=None,
        tick_base_minor=None,
        func_ax_lim=None,
        Di_ax_lim=None,
        thick=None,
        marker_size=None,
        func_marker_size=None,
        Di_marker_size=None,
        ref_m_picle=None,
        alpha=None,
        proj=None,
        rgba_bg=None,
        text_colour=None,
        text_ul=None,
        text_uc=None,
        text_ur=None,
        text_ll=None,
        text_lc=None,
        text_lr=None,
        Di_init_orbit=None,
        A1_po_inset=None,
        A1_inset_loc=None,
        fof_id_max=None,
        fof_m_min=None,
        fof_n_min=None,
        R_Roche=None,
        R_Hill=None,
        xz_size=None,
        zy_size=None,
        thick_xz=None,
        Fp_suffix=None,
        func_snapshot_add_misc=None,
        func_particle_params_extra=None,
        func_param_subsets_extra=None,
        func_reb_orbit_evol_extra=None,
    ):
        self.dir_save = dir_save
        self.category = category
        self.Di_mat_colour = Di_mat_colour
        self.Di_mat_highlt = Di_mat_highlt
        self.Di_param_unit_label = Di_param_unit_label
        self.unit_label_r_def = unit_label_r_def
        self.unit_label_m_def = unit_label_m_def
        self.Di_param_log = Di_param_log
        self.Di_param_lim = Di_param_lim
        self.Di_misc_colour = Di_misc_colour
        self.dpi = dpi
        self.fontsize = fontsize
        self.fontsize_text = fontsize_text
        self.cmap = cmap
        self.func_cmap_param = func_cmap_param
        self.Di_cmap_param = Di_cmap_param
        self.A1_ax_lim = A1_ax_lim
        self.tick_base = tick_base
        self.tick_base_minor = tick_base_minor
        self.func_ax_lim = func_ax_lim
        self.Di_ax_lim = Di_ax_lim
        self.thick = thick
        self.marker_size = marker_size
        self.func_marker_size = func_marker_size
        self.Di_marker_size = Di_marker_size
        self.ref_m_picle = ref_m_picle
        self.alpha = alpha
        self.proj = proj
        self.rgba_bg = rgba_bg
        self.text_colour = text_colour
        self.text_ul = text_ul
        self.text_uc = text_uc
        self.text_ur = text_ur
        self.text_ll = text_ll
        self.text_lc = text_lc
        self.text_lr = text_lr
        self.Di_init_orbit = Di_init_orbit
        self.A1_po_inset = A1_po_inset
        self.A1_inset_loc = A1_inset_loc
        self.fof_id_max = fof_id_max
        self.fof_m_min = fof_m_min
        self.fof_n_min = fof_n_min
        self.R_Roche = R_Roche
        self.R_Hill = R_Hill
        self.xz_size = xz_size
        self.zy_size = zy_size
        self.thick_xz = thick_xz
        self.Fp_suffix = Fp_suffix
        self.func_snapshot_add_misc = func_snapshot_add_misc
        self.func_particle_params_extra = func_particle_params_extra
        self.func_param_subsets_extra = func_param_subsets_extra
        self.func_reb_orbit_evol_extra = func_reb_orbit_evol_extra

        # Avoid recursive copying if this object is for an inset panel
        if self.A1_inset_loc is not None:
            A1_exclude = ["A1_po_inset"]
        else:
            A1_exclude = []

        # Copy
        if copy is not None:
            copy_object(self, copy, A1_exclude=A1_exclude)

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        # Set any unset dictionary entries to the defaults
        self.Di_mat_colour = copy_dict(self.Di_mat_colour, ut.Di_mat_colour)
        self.Di_mat_highlt = copy_dict(
            self.Di_mat_highlt, {mat: "#7711dd" for mat in self.Di_mat_colour.keys()}
        )
        if self.unit_label_r_def is None:
            self.unit_label_r_def = ParamUnitLabel(
                "r", "Distance", Ea.R, "$R_\oplus$", "R_E", " * Ea.R"
            )
        if self.unit_label_m_def is None:
            self.unit_label_m_def = ParamUnitLabel(
                "m", "Mass", Ea.M, "$M_\oplus$", "M_E", " * Ea.M"
            )
        Di_param_unit_label_def_custom = copy_dict(
            {
                # Set some non-SI defaults
                "m": ParamUnitLabel("m", copy=self.unit_label_m_def),
                "r": ParamUnitLabel("r", copy=self.unit_label_r_def),
                "x": ParamUnitLabel("x", copy=self.unit_label_r_def),
                "y": ParamUnitLabel("y", copy=self.unit_label_r_def),
                "z": ParamUnitLabel("z", copy=self.unit_label_r_def),
                "h": ParamUnitLabel("h", copy=self.unit_label_r_def),
                "a": ParamUnitLabel("a", copy=self.unit_label_r_def),
                "q": ParamUnitLabel("q", copy=self.unit_label_r_def),
                "Q": ParamUnitLabel("Q", copy=self.unit_label_r_def),
                "a_eq": ParamUnitLabel("a_eq", copy=self.unit_label_r_def),
            },
            Di_param_unit_label_def,
        )
        self.Di_param_unit_label = copy_dict(
            self.Di_param_unit_label, Di_param_unit_label_def_custom
        )
        Di_param_log_def = {
            "m": True,
            "u": True,
            "rho": True,
            "P": True,
        }
        self.Di_param_log = copy_dict(self.Di_param_log, Di_param_log_def)

        # Defaults
        if self.func_cmap_param is None:
            self.func_cmap_param = ut.func_cmap_param_def
        if self.dpi is None:
            self.dpi = 300
        if self.fontsize is None:
            self.fontsize = 30
        if self.fontsize_text is None:
            self.fontsize_text = 26
        if self.cmap is None:
            self.cmap = cmap_rbow
        if self.Di_cmap_param is None:
            self.Di_cmap_param = {}
        if self.Di_ax_lim is None:
            self.Di_ax_lim = {}
        if self.Di_marker_size is None:
            self.Di_marker_size = {}
        if self.alpha is None:
            self.alpha = 0.5
        if self.proj is None:
            self.proj = "xy"
        if self.rgba_bg is None:
            self.rgba_bg = (1.0, 1.0, 1.0, 0.85)
        if self.text_colour is None:
            self.text_colour = "k"
        if self.fof_n_min is None:
            self.fof_n_min = 50

        # Duplicate dictionaries with material IDs instead of names as keys
        if self.Di_mat_colour is not None:
            self.Di_id_colour = {}
            for mat, colour in self.Di_mat_colour.items():
                if mat[-2:] == "_2":
                    mat_id = woma.Di_mat_id[mat[:-2]] + ut.id_body
                else:
                    mat_id = woma.Di_mat_id[mat]
                self.Di_id_colour[mat_id] = colour
        if self.Di_mat_highlt is not None:
            self.Di_id_highlt = {}
            for mat, colour in self.Di_mat_highlt.items():
                if mat[-2:] == "_2":
                    mat_id = woma.Di_mat_id[mat[:-2]] + ut.id_body
                else:
                    mat_id = woma.Di_mat_id[mat]
                self.Di_id_highlt[mat_id] = colour

    def colours_from_data(self, param, A1_param, Di_cmap=None):
        """Convert parameter data to colours for plotting.

        Parameters
        ----------
        param : str
            The parameter used to set the colours, e.g. a particle property.

        A1_param : [?]
            The parameter array (plotting units).

        Di_cmap : {} (opt.)
            A dictionary with colour-map info required for some parameters.
            "cmap" : Colormap
                A colour map.

            "norm" : Normalize
                A normalisation object.

            #
            # Special cases can require extra entries, e.g. a second colour map.
            #

        Returns
        -------
        A1_rgba : [[float]]
            The array of RGBA colours.
        """
        if param == "mat_id":
            # Initialise colour array
            A1_rgba = np.empty((len(A1_param), 4), dtype=object)

            # Set each material
            for mat_id, colour in self.Di_id_colour.items():
                A1_sel_mat = np.where(A1_param == mat_id)[0]

                if len(A1_sel_mat) > 0:
                    A1_rgba[A1_sel_mat] = list(mpl.colors.to_rgba(colour))

            if None in A1_rgba:
                raise Exception(
                    "Not all materials converted to colours ", np.unique(A1_param)
                )
        elif param == "fof_id":
            # Initialise colour array
            A1_rgba = np.empty((len(A1_param), 4), dtype=object)

            # Select large groups, individual colours
            A1_fof_id_plot = np.unique(A1_param)
            A1_fof_id_plot = A1_fof_id_plot[A1_fof_id_plot <= self.fof_id_max]
            A1_c_fof = self.cmap(np.linspace(0, 1, len(A1_fof_id_plot)))

            # Colour each selected large group
            for fof_id, colour in zip(A1_fof_id_plot, A1_c_fof):
                A1_sel_fof = np.where(A1_param == fof_id)[0]
                A1_rgba[A1_sel_fof] = list(mpl.colors.to_rgba(colour))

            # Colour all smaller groups
            A1_sel_small = np.where(
                (self.fof_id_max < A1_param) & (A1_param != ut.fof_id_none)
            )[0]
            A1_rgba[A1_sel_small] = list(
                mpl.colors.to_rgba(self.Di_misc_colour["fof_small"])
            )

            # Colour all no-group particles
            A1_sel_none = np.where(A1_param == ut.fof_id_none)[0]
            A1_rgba[A1_sel_none] = list(
                mpl.colors.to_rgba(self.Di_misc_colour["fof_none"])
            )
        elif param == "Q_v_inf":
            cmap = Di_cmap["cmap"]
            norm = Di_cmap["norm"]
            A1_sel = Di_cmap["A1_sel"]
            cmap_2 = Di_cmap["cmap_2"]
            norm_2 = Di_cmap["norm_2"]
            A1_sel_2 = Di_cmap["A1_sel_2"]

            # Initialise colour array, then set each subset
            A1_rgba = np.empty((len(A1_param), 4))
            A1_rgba[A1_sel] = cmap(norm(A1_param[A1_sel]))
            A1_rgba[A1_sel_2] = cmap_2(norm_2(A1_param[A1_sel_2]))
        else:
            cmap = Di_cmap["cmap"]
            norm = Di_cmap["norm"]

            # RGBA colours from colour map
            A1_rgba = cmap(norm(A1_param))

        # Opacity
        A1_rgba[:, 3] = self.alpha

        # Ensure no invalid values
        A1_rgba[A1_rgba < 0] = 0
        A1_rgba[A1_rgba > 1] = 1

        return A1_rgba

    def mark_inset_panel(self, ax, ax_in):
        """Mark an inset zoom in or out panel.

        Parameters
        ----------
        ax : Axes
            The main figure axes.

        ax_in : Axes
            The inset panel axes.
        """
        # Is the inset a zoom-in or a zoom-out?
        x_min, x_max = ax.get_xlim()
        x_min_in, x_max_in = ax_in.get_xlim()
        is_zoom_in = x_min < x_min_in and x_max > x_max_in

        # Zoom in: mark the inset on the main axes, and connect two corners
        if is_zoom_in:
            # Corners to connect: 1=top right, 2=top left, 3=bottom left, 4=bottom right
            is_left = self.A1_inset_loc[0] < 0.1
            is_top = self.A1_inset_loc[1] > 0.1
            if is_left * is_top:
                # Top left or bottom right
                loc_1 = 1
                loc_2 = 3
            else:
                # Top right or bottom left
                loc_1 = 2
                loc_2 = 4

            mark_inset(
                ax,
                ax_in,
                loc1=loc_1,
                loc2=loc_2,
                fc="none",
                ec="k",
                lw=0.5,
                alpha=0.7,
            )

        # Zoom out: mark the main axes on the inset
        else:
            # Main axis coordinates with minor padding
            d_pad = 0.03
            y_min, y_max = ax.get_ylim()
            dx = d_pad * (x_max - x_min)
            x_min -= dx
            x_max += dx
            dy = d_pad * (y_max - y_min)
            y_min -= dy
            y_max += dy

            ax_in.plot(
                [x_min, x_max, x_max, x_min, x_min],
                [y_min, y_min, y_max, y_max, y_min],
                c="k",
                lw=0.8,
                alpha=0.6,
            )

    def plot_init_orbit(self, o_0, ax=None):
        """Plot an initial orbit.

        Parameters
        ----------
        o_0 : Orbit
            The initial orbit.

        ax : Axes (opt.)
            The axes to plot on.

        Required self.Di_init_orbit
        ---------------------------
        "colour", "lw", "alpha" : str, float, float
            The colour, linewidth and opacity for the orbit.

        "ls", or "ls_step_on", "ls_step_off" : float, float
            The linestyle parameters for the orbit.
        """
        if ax is None:
            ax = plt.gca()

        # Placeholders
        if "ls" not in self.Di_init_orbit.keys():
            self.Di_init_orbit["ls"] = None
        if "ls_step_on" not in self.Di_init_orbit.keys():
            self.Di_init_orbit["ls_step_on"] = None
            self.Di_init_orbit["ls_step_off"] = None

        plot_orbit(
            o_0,
            A1_pos_0=[0, 0, 0],
            unit=self.r_unit,
            proj=self.proj,
            c=self.Di_init_orbit["colour"],
            lw=self.Di_init_orbit["lw"],
            alpha=self.Di_init_orbit["alpha"],
            ls=self.Di_init_orbit["ls"],
            ls_step_on=self.Di_init_orbit["ls_step_on"],
            ls_step_off=self.Di_init_orbit["ls_step_off"],
            A1_ax_lim="auto",
            zorder=-999,
            ax=ax,
        )

    @property
    def r_unit(self):
        return self.Di_param_unit_label["r"].unit

    @property
    def r_unit_label(self):
        return self.Di_param_unit_label["r"].unit_label


class InitProf(woma.Planet):
    """Planet profiles for initial conditions.

    Based on the WoMa Planet class, with extra bits for e.g. project management,
    setting arguments for woma functions, and storing output data for reference.

    Parameters
    ----------
    name : str
        The planet name. Nominally prefixed with "prof_".

    copy : str
        The name of a different InitProf object from which to copy all
        attributes that are not provided as other arguments here.

    dir_proj : str
        The file path to the base directory for the project. Expected to contain
        a sub directory called init_cond/.

    Fp_save : str
        The file path for the saved profiles. Defaults to
        dir_proj/init_cond/name.hdf5.

    po : PlotOptions
        An object with settings for plotting functions.

    planet : InitProf
        For spinning planets, the base spherical planet profile.

    period : float
        For spinning planets, the period (h).

    #
    # woma.Planet parameters...
    #

    #
    # woma.SpinPlanet parameters...
    #

    R_min, R_max : float
        The minimum and maximum radii to try (m).

    tol : float
        The tolerance for finding unknown parameters as a fractional difference
        between two consecutive iterations.

    tol_M_tweak : float
        The tolerance for tweaking the mass to avoid density peaks at the
        centre; the relative difference between consecutive masses.

    tol_M_layer : float
        The tolerance for the desired outer-layer mass when tweaking the surface
        pressure to integrate outwards.

    num_attempt : int
        The maximum number of iteration attempts.

    num_prof_spin : int
        For spinning planets, the number of grid points used in the 1D
        equatorial and polar profiles, i.e. the number of nested spheroids used
        to model the spinning planet.

    rho_min, P_min : float
        The minimum density and/or pressure (must be >= 0) at which the new
        layer will stop.

    L_max : float
        A reference maximum spin angular momentum, e.g. as a plotting unit.

    Di_po_edit : {}
        A dictionary of parameter names and values to edit the PlotOptions
        object.
    """

    # Dictionary of all instances: {self.name : self}
    _Di_init_prof = {}

    def __init__(
        self,
        name,
        copy=None,
        dir_proj=None,
        Fp_save=None,
        po=None,
        planet=None,
        period=None,
        # woma.Planet output parameters
        A1_mat_layer=None,
        A1_T_rho_type=None,
        P_s=None,
        T_s=None,
        rho_s=None,
        M=None,
        R=None,
        A1_M_layer=None,
        A1_R_layer=None,
        A1_Z_layer=None,
        A1_idx_layer=None,
        P_0=None,
        T_0=None,
        rho_0=None,
        P_1=None,
        T_1=None,
        rho_1=None,
        P_2=None,
        T_2=None,
        rho_2=None,
        I_MR2=None,
        L=None,
        num_prof=None,
        verbosity=None,
        # woma.SpinPlanet output parameters
        num_prof_spin=None,
        f_iter=None,
        tol_density_profile=None,
        tol_layer_masses=None,
        num_attempt_1=None,
        num_attempt_2=None,
        # Additional parameters
        R_min=None,
        R_max=None,
        tol=None,
        tol_M_tweak=None,
        tol_M_layer=None,
        num_attempt=None,
        rho_min=None,
        P_min=None,
        L_max=None,
        Di_po_edit=None,
    ):
        # Set initial attributes before calling the main woma.Planet's init, to
        # allow copying of unset attributes first
        self.name = name
        self.copy = copy
        self.dir_proj = dir_proj
        self.Fp_save = Fp_save
        self.po = po
        self.planet = planet
        self.period = period
        # woma.Planet parameters
        self.A1_mat_layer = A1_mat_layer
        self.A1_T_rho_type = A1_T_rho_type
        self.P_s = P_s
        self.T_s = T_s
        self.rho_s = rho_s
        self.M = M
        self.R = R
        self.A1_M_layer = A1_M_layer
        self.A1_R_layer = A1_R_layer
        self.A1_Z_layer = A1_Z_layer
        self.A1_idx_layer = A1_idx_layer
        self.P_0 = P_0
        self.T_0 = T_0
        self.rho_0 = rho_0
        self.P_1 = P_1
        self.T_1 = T_1
        self.rho_1 = rho_1
        self.P_2 = P_2
        self.T_2 = T_2
        self.rho_2 = rho_2
        self.I_MR2 = I_MR2
        self.L = L
        self.num_prof = num_prof
        self.verbosity = verbosity
        # woma.SpinPlanet parameters
        self.num_prof_spin = num_prof_spin
        self.f_iter = f_iter
        self.tol_density_profile = tol_density_profile
        self.tol_layer_masses = tol_layer_masses
        self.num_attempt_1 = num_attempt_1
        self.num_attempt_2 = num_attempt_2
        # Additional parameters
        self.R_min = R_min
        self.R_max = R_max
        self.tol = tol
        self.tol_M_tweak = tol_M_tweak
        self.tol_M_layer = tol_M_layer
        self.num_attempt = num_attempt
        self.rho_min = rho_min
        self.P_min = P_min
        self.L_max = L_max
        self.Di_po_edit = Di_po_edit

        self.init_done = False

        # Add to the dictionary
        if self.name is not None:
            InitProf._Di_init_prof[self.name] = self

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        # Copy
        if self.copy is not None:
            copy = InitProf._Di_init_prof[self.copy]
            copy.init()
            copy_object(
                self,
                copy,
                A1_exclude=["Fp_save", "Di_po_edit"],
                A1_set_eq=["po", "planet"],
            )

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        # Defaults
        if self.Fp_save is None:
            self.Fp_save = "%s/init_cond/%s.hdf5" % (self.dir_proj, self.name)
        if self.po is None and self.period is None:
            self.po = PlotOptions()
        if self.num_prof is None:
            self.num_prof = 10000
        if self.verbosity is None:
            self.verbosity = 1
        if self.num_prof_spin is None:
            self.num_prof_spin = 3000
        if self.f_iter is None:
            self.f_iter = 0.002
        if self.tol_density_profile is None:
            self.tol_density_profile = 0.02
        if self.tol_layer_masses is None:
            self.tol_layer_masses = 0.01
        if self.num_attempt_1 is None:
            self.num_attempt_1 = 15
        if self.num_attempt_2 is None:
            self.num_attempt_2 = 5
        if self.tol is None:
            self.tol = 0.0001
        if self.tol_M_tweak is None:
            self.tol_M_tweak = 1e-7
        if self.tol_M_layer is None:
            self.tol_M_layer = 0.005
        if self.num_attempt is None:
            self.num_attempt = 40
        if self.rho_min is None:
            self.rho_min = 0
        if self.P_min is None:
            self.P_min = 0

        # Main init
        if self.period is None:
            woma.Planet.__init__(
                self,
                name=self.name,
                A1_mat_layer=self.A1_mat_layer,
                A1_T_rho_type=self.A1_T_rho_type,
                P_s=self.P_s,
                T_s=self.T_s,
                rho_s=self.rho_s,
                M=self.M,
                R=self.R,
                A1_M_layer=self.A1_M_layer,
                A1_R_layer=self.A1_R_layer,
                A1_idx_layer=self.A1_idx_layer,
                P_0=self.P_0,
                T_0=self.T_0,
                rho_0=self.rho_0,
                P_1=self.P_1,
                T_1=self.T_1,
                rho_1=self.rho_1,
                P_2=self.P_2,
                T_2=self.T_2,
                rho_2=self.rho_2,
                I_MR2=self.I_MR2,
                num_prof=self.num_prof,
                load_file=None,
                verbosity=self.verbosity,
            )
        else:
            self.planet.init()

            # Copy unset parameters from the base spherical planet
            if self.dir_proj is None:
                self.dir_proj = self.planet.dir_proj
            if self.po is None:
                self.po = self.planet.po
            if not hasattr(self, "A1_mat_id_layer") or self.A1_mat_id_layer is None:
                self.A1_mat_id_layer = self.planet.A1_mat_id_layer
            if not hasattr(self, "A1_mat_layer") or self.A1_mat_layer is None:
                self.A1_mat_layer = self.planet.A1_mat_layer

        # Overrides
        if self.Di_po_edit is not None:
            # Make a clean copy to edit
            self.po = deepcopy(self.po)

            for var in self.Di_po_edit.keys():
                # self.po.var = self.Di_po_edit["var"]
                exec("self.po." + var + ' = self.Di_po_edit["' + var + '"]')

    def load_profiles(self):
        """Load the profile data (see WoMa)."""
        if self.period is None:
            self.load(self.Fp_save)
        else:
            print('Loading "%s"... ' % self.Fp_save[-58:], end="", flush=True)
            with h5py.File(self.Fp_save, "r") as f:
                (
                    self.A1_R,
                    self.A1_Z,
                    self.A1_m,
                    self.A1_rho,
                    self.A1_T,
                    self.A1_P,
                    self.A1_u,
                    self.A1_mat_id,
                ) = woma.misc.io.multi_get_spin_planet_data(
                    f, ["R", "Z", "m", "rho", "T", "P", "u", "mat_id"]
                )
            print("Done")


class InitCond:
    """Initial conditions including placing particles (for a single body).

    Parameters
    ----------
    name : str
        The planet name. Nominally prefixed with "init_".

    copy : str
        The name of a different InitCond object from which to copy all
        attributes that are not provided as other arguments here.

    ip : InitProf
        The planet profile object.

    num_picle_des : int
        The desired number of particles to place.

    Fp_save : str
        The file path for the saved initial conditions. Defaults to
        ip.dir_proj/name.hdf5.

    N_ngb : int
        The number of neighbours used to estimate the SPH smoothing lengths and
        densities.

    file_to_SI : Conversions
        The unit conversion object for the file units into SI. See jkeger.py.

    boxsize : float
        The extent of the simulation box (m).

    verbosity : int
        The verbosity to control printed output.

    num_picle : int
        The resulting number of particles. May differ slightly from num_picle_des.

    M : float
        The total mass of particles (kg).

    R_s, Z_s, R_s_sph : float
        The outer "surface" radius (if a spinning planet, then the equatorial
        and polar radii, and _sph the spherical base profile's radius), e.g.
        ignoring any atmosphere (m).

    min_sep : float
        The minimum inter-particle separation (m).

    Attributes
    ----------
    po : PlotOptions
        Settings for plotting functions, inherited from ip.

    SI_to_file : Conversions
        Inverse of file_to_SI.

    m_picle : float
        The mean particle mass.

    v_esc : float
        The escape speed from the surface (m s^-1).
    """

    # Dictionary of all instances: {self.name : self}
    _Di_init_cond = {}

    def __init__(
        self,
        name,
        copy=None,
        ip=None,
        num_picle_des=None,
        Fp_save=None,
        N_ngb=None,
        file_to_SI=None,
        boxsize=None,
        verbosity=None,
        num_picle=None,
        M=None,
        R_s=None,
        Z_s=None,
        R_s_sph=None,
        min_sep=None,
    ):
        # Set instead of using woma's init that places the particles immediately
        self.name = name
        self.copy = copy
        self.ip = ip
        self.num_picle_des = num_picle_des
        self.Fp_save = Fp_save
        self.N_ngb = N_ngb
        self.file_to_SI = file_to_SI
        self.boxsize = boxsize
        self.verbosity = verbosity
        self.num_picle = num_picle
        self.M = M
        self.R_s = R_s
        self.Z_s = Z_s
        self.R_s_sph = R_s_sph
        self.min_sep = min_sep

        self.init_done = False

        # Add to the dictionary
        if self.name is not None:
            InitCond._Di_init_cond[self.name] = self

    def init(self):
        # Don't init twice
        if self.init_done:
            return
        self.init_done = True

        # Copy
        if self.copy is not None:
            copy = InitCond._Di_init_cond[self.copy]
            copy.init()
            copy_object(
                self,
                copy,
                A1_exclude=[
                    "Fp_save",
                    "num_picle",
                    "M",
                    "R_s",
                    "Z_s",
                    "R_s_sph",
                    "min_sep",
                ],
                A1_set_eq=["ip"],
            )

        # Overwrite "None" to None
        set_selected_none_attributes(self)

        self.ip.init()

        # Defaults
        if self.Fp_save is None:
            self.Fp_save = "%s/init_cond/%s.hdf5" % (self.ip.dir_proj, self.name)
        if self.N_ngb is None:
            self.N_ngb = 48
        if self.verbosity is None:
            self.verbosity = 1
        if self.R_s is None:
            # Set the surface from the profile, excluding a H--He atmosphere
            if self.ip.A1_mat_layer[-1] in ["HM80_HHe"]:
                self.R_s = self.ip.A1_R_layer[-2]
                if self.ip.period is not None:
                    if self.Z_s is None:
                        self.Z_s = self.ip.A1_Z_layer[-2]
                    if self.R_s_sph is None:
                        self.R_s_sph = self.ip.planet.A1_R_layer[-2]
            else:
                self.R_s = self.ip.A1_R_layer[-1]
                if self.ip.period is not None:
                    if self.Z_s is None:
                        self.Z_s = self.ip.A1_Z_layer[-1]
                    if self.R_s_sph is None:
                        self.R_s_sph = self.ip.planet.A1_R_layer[-1]

        # Derived attributes
        self.po = self.ip.po
        self.SI_to_file = self.file_to_SI.inv()

    @property
    def m_picle(self):
        return self.M / self.num_picle

    @property
    def v_esc(self):
        return v_esc_from_M_R(self.M, self.R_s)


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR classes.py  ====\n")
