"""
All my useful python things. Standard constants and generic utilities etc. that
aren't specific to individual projects.

Note that some blocks are indented in `if True:` statements purely for
convenient block-folding in the code editor.
"""
import numpy as np
import sys
import os
from numba import njit
import pickle
import random
import tempfile
import inspect
import warnings
import fileinput
import shutil
import traceback
from copy import deepcopy
from time import sleep, mktime
import datetime
import colorsys
import string
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import sem, binned_statistic, binned_statistic_2d
from scipy.optimize import curve_fit, fsolve, minimize_scalar
from scipy.interpolate import make_interp_spline, BSpline, CubicSpline
from scipy.spatial.transform import Rotation

# Path to this file's directory
path = os.path.dirname(os.path.realpath(__file__))

# Create a temporary dir for matplotlib lock files to avoid conflicts
tmp_dir_mpl = tempfile.mkdtemp()
os.environ["MPLCONFIGDIR"] = tmp_dir_mpl
import matplotlib as mpl

# Use a non-interactive backend for saving figures
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse, Arc, FancyArrowPatch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

try:
    # Import WoMa if available
    sys.path.append(os.path.join(path, "../WoMa/"))
    import woma
except ImportError:
    pass

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

from orbits import Orbit, OrbitSet


# ========
# Constants and conversions
# ========
pi = np.pi
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
R_gas = 8.3145  # Gas constant (J K^-1 mol^-1)
k_B = 1.38065e-23  # Boltzmann constant (kg m^2 s^-2 K^-1)

# Unit conversions
class Conversions:
    """Class to store conversions from one set of units to another, derived
    using the base mass-, length-, and time-unit relations.

    Usage e.g.
    ----------
    cgs_to_SI = Conversions(m=1e-3, l=1e-2, t=1)
    SI_to_cgs = cgs_to_SI.inv()

    rho_SI = rho_cgs * cgs_to_SI.rho
    G_cgs = 6.67e-11 * SI_to_cgs.G

    Parameters
    ----------
    m : float
        Value to convert mass from the first units to the second.

    l : float
        Value to convert length from the first units to the second.

    t : float
        Value to convert time from the first units to the second.

    Attributes (all floats)
    ----------
    m           Mass
    l           Length
    t           Time
    v           Velocity
    a           Acceleration
    rho         Density
    drho_dt     Rate of change of density
    P           Pressure
    u           Specific energy
    du_dt       Rate of change of specific energy
    E           Energy
    s           Specific entropy
    G           Gravitational constant
    L           Angular momentum
    """

    def __init__(self, m, l, t):
        # Input conversions
        self.m = m
        self.l = l
        self.t = t
        # Derived conversions
        self.v = l * t**-1
        self.a = l * t**-2
        self.rho = m * l**-3
        self.drho_dt = m * l**-4
        self.P = m * l**-1 * t**-2
        self.u = l**2 * t**-2
        self.du_dt = l**2 * t**-3
        self.E = m * l**2 * t**-2
        self.s = l**2 * t**-2
        self.G = m**-1 * l**3 * t**-2
        self.L = m * l**2 * t**-1

    def inv(self):
        """Return the inverse to this conversion"""
        return Conversions(1 / self.m, 1 / self.l, 1 / self.t)


# Typical units
if True:
    SI_to_SI = Conversions(m=1, l=1, t=1)  # No-op

    # Standard
    cgs_to_SI = Conversions(m=1e-3, l=1e-2, t=1)
    SI_to_cgs = cgs_to_SI.inv()
    # Angles
    deg_to_rad = pi / 180
    rad_to_deg = 1 / deg_to_rad
    deg_to_arcmin = 60
    arcmin_to_deg = 1 / deg_to_arcmin
    deg_to_arcsec = 60**2
    arcsec_to_deg = 1 / deg_to_arcsec
    rad_to_arcmin = rad_to_deg * deg_to_arcmin
    arcmin_to_rad = 1 / rad_to_arcmin
    rad_to_arcsec = rad_to_deg * deg_to_arcsec
    arcsec_to_rad = 1 / rad_to_arcsec
    # Time
    hour_to_s = 60**2  # = 3.60e3 s
    s_to_hour = 1 / hour_to_s  # = 2.77e-4 h
    day_to_s = 24 * hour_to_s  # = 8.64e4 s
    s_to_day = 1 / day_to_s  # = 1.16e-5 days
    week_to_s = 7 * day_to_s  # = 6.05e5 s
    s_to_week = 1 / week_to_s  # = 1.65e-6 weeks
    yr_to_s = 365.25 * day_to_s  # = 3.15e7 s
    s_to_yr = 1 / yr_to_s  # = 3.17e-8 yr
    sidereal_yr_to_s = 31558149.7635
    s_to_sidereal_yr = 1 / sidereal_yr_to_s
    # Distance
    mile_to_m = 1609.34
    m_to_mile = 1 / mile_to_m  # = 6.21e-4 miles
    au_to_m = 1.4959787e11
    m_to_au = 1 / au_to_m  # = 6.68e-12 AU
    pc_to_m = 3.085677581e16
    m_to_pc = 1 / pc_to_m
    # Speed
    mph_to_mps = mile_to_m / hour_to_s  # = 0.447 m/s
    mps_to_mph = 1 / mph_to_mps  # = 2.237 mph
    kmph_to_mps = 1e3 / hour_to_s  # = 0.278 m/s
    mps_to_kmph = 1 / kmph_to_mps  # = 3.6 km/h
    # Pressure
    bar_to_Pa = 1e5
    Pa_to_bar = 1 / bar_to_Pa  # = 1e-5 bar
    Mbar_to_Pa = bar_to_Pa * 1e6  # = 1e11 Pa
    Pa_to_Mbar = 1 / Mbar_to_Pa  # = 1e-11 Mbar
    bar_to_Ba = bar_to_Pa * SI_to_cgs.P  # = 1e6 Ba
    Ba_to_bar = 1 / bar_to_Ba  # = 1e-6 bar
    Mbar_to_Ba = Mbar_to_Pa * SI_to_cgs.P  # = 1e12 Ba
    Ba_to_Mbar = 1 / Mbar_to_Ba  # = 1e-12 Mbar
    # Mass
    amu_to_kg = 1.6605e-27
    kg_to_amu = 1 / amu_to_kg  # = 6.02e26 amu

    # ========
    # Misc
    # ========
    # n-sigma amounts and percentiles
    sigma_1 = 0.6827
    sigma_2 = 0.9545
    sigma_3 = 0.9974
    sigma_1_low = 0.5 * (1 - sigma_1)
    sigma_1_high = 0.5 * (1 + sigma_1)
    sigma_2_low = 0.5 * (1 - sigma_2)
    sigma_2_high = 0.5 * (1 + sigma_2)
    sigma_3_low = 0.5 * (1 - sigma_3)
    sigma_3_high = 0.5 * (1 + sigma_3)

    # Earth--Moon angular momentum
    L_EM = 3.5e34  # kg m^2 s^-1


# ========
# Parameter labels and units etc
# ========
Di_param_label = {
    "id": r"Particle ID",
    "m": r"Mass",
    "mat_id": r"Material ID",
    "fof_id": r"FoF group ID",
    "r": r"Radial distance",
    "x": r"$x$",
    "y": r"$y$",
    "z": r"$z$",
    "v": r"Speed",
    "v_x": r"$v_x$",
    "v_y": r"$v_y$",
    "v_z": r"$v_z$",
    "rho": r"Density",
    "h": r"Smoothing length",
    "u": r"Specific internal energy",
    "P": r"Pressure",
    "T": r"Temperature",
    "s": r"Specific entropy",
    "phi": r"Gravitational potential",
    "t": r"Time",
    "E": r"Energy",
    "L": r"Angular momentum",
    "b": r"Impact parameter",
    "beta": r"Impact angle",
    "a": r"Semi-major axis",
    "e": r"Eccentricity",
    "i": r"Inclination",
    "q": r"Periapsis",
    "Q": r"Apoapsis",
    "Omega": r"Longitude of ascending node",
    "pomega": r"Longitude of periapsis",
    "omega": r"Argument of periapsis",
    "period": r"Period",
    "a_eq": r"Equivalent circular orbit",
}
Di_param_symbol = {
    "id": r"ID",
    "m": r"$m$",
    "r": r"$r$",
    "x": r"$x$",
    "y": r"$y$",
    "z": r"$z$",
    "v": r"$v$",
    "v_x": r"$v_x$",
    "v_y": r"$v_y$",
    "v_z": r"$v_z$",
    "rho": r"$\rho$",
    "h": r"$h$",
    "u": r"$u$",
    "P": r"$P$",
    "T": r"$T$",
    "s": r"$s$",
    "phi": r"$\phi$",
    "t": r"$t$",
    "E": r"$E$",
    "L": r"$L$",
    "b": r"$b$",
    "beta": r"$\beta$",
    "a": r"$a$",
    "e": r"$e$",
    "i": r"$i$",
    "q": r"$q$",
    "Q": r"$Q$",
    "Omega": r"$\Omega$",
    "pomega": r"$\varpi$",
    "omega": r"$\omega$",
    "period": r"$T$",
    "a_eq": r"$a_{\rm eq}$",
}


class ParamUnitLabel:
    """Class for simple parameter units and labels, etc.

    Parameters
    ----------
    param : str
        The parameter name, e.g. "m".

    label : str (opt.)
        The main label text (defaults from Di_param_label, required otherwise),
        e.g. "Mass".

    unit : float (opt.)
        The plotting-etc units (default 1 for SI), e.g. Ea.M.

    unit_label : str (opt.)
        The units for latex plotting etc, e.g. r"$M_\oplus$".

    unit_print : str (opt.)
        The units for screen printing, e.g. "M_E".

    unit_code : str (opt.)
        The units as an executable inline code string, e.g. "* Ea.M".

    symbol : str (opt.)
        The parameter symbol for latex plotting etc, e.g. r"$m$".

    copy : ParamUnitLabel (opt.)
        If provided, then copy any not-provided inputs from this object.
    """

    def __init__(
        self,
        param,
        label=None,
        unit=None,
        unit_label=None,
        unit_print=None,
        unit_code=None,
        symbol=None,
        copy=None,
    ):
        self.param = param
        self.label = label
        self.unit = unit
        self.unit_label = unit_label
        self.unit_print = unit_print
        self.unit_code = unit_code
        self.symbol = symbol
        self.copy = copy

        # Copy
        if self.copy is not None:
            copy_object(self, self.copy, A1_exclude=["param", "label"])

        # Defaults
        if self.label is None:
            # Required if not available as a standard default
            assert self.param in Di_param_label.keys()
            self.label = Di_param_label[self.param]
        if self.unit == None:
            self.unit = 1
        if self.unit == 1:
            self.unit_code = ""
        if self.symbol is None and self.param in Di_param_symbol.keys():
            self.symbol = Di_param_symbol[self.param]

        # ========
        # Preset defaults
        # ========
        # Mass
        if self.param == "m":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "kg"

        # Distance
        elif self.param in ["r", "x", "y", "z", "v", "h", "a", "q", "Q", "a_eq"]:
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "m"

        # Speed
        elif self.param in ["v", "v_x", "v_y", "v_z"]:
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = r"m s$^{-1}$"
                if self.unit_print is None:
                    self.unit_print = r"m s^-1"

        # Energy per mass
        elif self.param in ["u", "phi"]:
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = r"$\rm{J~kg}^{-1}$"
                if self.unit_print is None:
                    self.unit_print = "J kg^-1"

        # Misc
        elif self.param == "rho":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = r"$\rm{kg~m}^{-3}$"
                if self.unit_print is None:
                    self.unit_print = "kg m^-3"
        elif self.param == "P":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "Pa"
        elif self.param == "T":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "K"
        elif self.param == "s":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = r"$\rm{J~kg}^{-1}\rm{~K}^{-1}$"
                if self.unit_print is None:
                    self.unit_print = "J kg^-1 K^-1"
        elif self.param == "t":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "s"
        elif self.param == "E":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "J"
        elif self.param == "L":
            if self.unit == 1:
                if self.unit_label is None:
                    self.unit_label = "$\rm{kg~m}^{2}\rm{~s}^{-1}$"
                if self.unit_print is None:
                    self.unit_print = "kg m^2 s^-1"

        # Degrees
        elif self.unit is deg_to_rad:
            if self.unit_label is None:
                self.unit_label = r"$^\circ$"
            if self.unit_print is None:
                self.unit_print = "deg"

        # ========
        # Leftovers
        # ========
        # No units
        if self.unit is None or self.unit_label is None:
            self.unit = 1
            self.unit_label = ""
            self.unit_print = ""

        # Default print
        if self.unit_print is None:
            self.unit_print = self.unit_label

    @property
    def label_unit(self):
        """A label with units for latex plotting etc, e.g. r"Mass $(M_\oplus)$"."""
        if self.unit_label == "":
            return self.label
        else:
            return r"%s (%s)" % (self.label, self.unit_label)


# Default labels and SI units for standard parameters
Di_param_unit_label = {param: ParamUnitLabel(param) for param in Di_param_label.keys()}

# Override default angle units to degrees
for param in ["i", "beta", "Omega", "pomega", "omega"]:
    Di_param_unit_label[param] = ParamUnitLabel(param, unit=deg_to_rad)


def roman(number, capital=False):
    """Convert an integer to roman numerals."""
    number = int(number)

    A1_int = [1, 4, 5, 9, 10, 40, 50, 90, 100, 400, 500, 900, 1000]
    if capital:
        A1_roman = [
            "I",
            "IV",
            "V",
            "IX",
            "X",
            "XL",
            "L",
            "XC",
            "C",
            "CD",
            "D",
            "CM",
            "M",
        ]
    else:
        A1_roman = [
            "i",
            "iv",
            "v",
            "ix",
            "x",
            "xl",
            "l",
            "xc",
            "c",
            "cd",
            "d",
            "cm",
            "m",
        ]

    i = len(A1_int) - 1
    roman = ""
    while number > 0:
        div = number // A1_int[i]
        number %= A1_int[i]

        while div > 0:
            roman += A1_roman[i]
            div -= 1
        i -= 1

    return roman


# Enumerating
A1_enum_1 = ["%d" % i for i in range(99)]
A1_enum_i = [roman(i) for i in range(99)]
A1_enum_I = [roman(i, capital=True) for i in range(99)]
A1_enum_a = list(string.ascii_lowercase)
A1_enum_A = list(string.ascii_uppercase)


# ========
# Avoid numpy precision errors
# ========
def mean(*args, **kwargs):
    return np.mean(*args, **kwargs, dtype=np.float64)


def sum(*args, **kwargs):
    return np.sum(*args, **kwargs, dtype=np.float64)


def cumsum(*args, **kwargs):
    return np.cumsum(*args, **kwargs, dtype=np.float64)


# ========
# Maths etc
# ========
def delta_ij(i, j):
    """Kroneckar delta."""
    if i == j:
        return 1
    else:
        return 0


def log_modulus_transform(x):
    """Return the log-modulus transformation of x, preserves 0 and -ve values"""
    return np.sign(x) * np.log10(abs(x) + 1)


def rotation_matrix(A1_u, theta):
    """Return the matrix to 3D rotate by an angle theta (rad) about a vector u."""
    A1_u = np.array(A1_u, dtype=float)
    A1_u /= np.linalg.norm(A1_u)

    cos_t = np.cos(theta)
    omcos_t = 1 - cos_t
    sin_t = np.sin(theta)

    return np.array(
        [
            [
                cos_t + A1_u[0] ** 2 * omcos_t,
                A1_u[0] * A1_u[1] * omcos_t - A1_u[2] * sin_t,
                A1_u[0] * A1_u[2] * omcos_t + A1_u[1] * sin_t,
            ],
            [
                A1_u[0] * A1_u[1] * omcos_t + A1_u[2] * sin_t,
                cos_t + A1_u[1] ** 2 * omcos_t,
                A1_u[1] * A1_u[2] * omcos_t - A1_u[0] * sin_t,
            ],
            [
                A1_u[0] * A1_u[2] * omcos_t - A1_u[1] * sin_t,
                A1_u[1] * A1_u[2] * omcos_t + A1_u[0] * sin_t,
                cos_t + A1_u[2] ** 2 * omcos_t,
            ],
        ]
    )


def find_index_and_interp(x, A1_x):
    """Return the index and interpolation factor of a value in an array.

    Allows x outside A1_x. If so then intp will be < 0 or > 1.

    Parameters
    ----------
    x : float
        The value to find.

    A1_x : [float]
        The array to search.

    Returns
    -------
    idx : int
        The index of the last array element smaller than the value.

        0               If x is below A1_x.
        len(A1_x) - 2   If x is above A1_x.

    intp : float
        The interpolation factor for how far the values is from the
        indexed array value to the next.

        < 0     If x is below A1_x.
        > 1     If x is above A1_x.
    """
    idx = np.searchsorted(A1_x, x) - 1

    # Return error values if outside the array
    if idx == -1:
        idx = 0
    elif idx >= len(A1_x) - 1:
        idx = len(A1_x) - 2

    intp = (x - A1_x[idx]) / (A1_x[idx + 1] - A1_x[idx])

    return idx, intp


@njit
def root_sum_sq(A1_x):
    """Return the root of the summed squares of all x in A1_x (incl. 1D x)."""
    ans = A1_x[0] ** 2
    for x in A1_x[1:]:
        ans += x**2

    return np.sqrt(ans)


def idx_closest(A1_x, x):
    """Return the index of the element of A1_x closest to x."""
    return (np.abs(A1_x - x)).argmin()


@njit
def mod_angle(a):
    """Set an angle between 0 and 2pi"""
    return np.mod(2 * pi + np.mod(a, 2 * pi), 2 * pi)


@njit
def arccos2(o, a, disambiguator=1):
    """Return arccos(o/a), using the disambiguator to set the quadrant."""
    cosine = o / a
    # In normal bounds
    if cosine > -1 and cosine < 1:
        angle = np.arccos(cosine)
        if disambiguator < 0:
            angle = -angle
    # Exceptions
    elif cosine <= -1:
        angle = pi
    else:
        angle = 0

    return angle


def round_to_nearest(value, round, prec=7):
    """Round the value to the nearest arbitrary amount."""
    return np.round((np.round(value / round)) * round, prec)


def inertia_tensor(A1_x, A1_y, A1_z, A1_m):
    """Return the inertia tensor for a set of points.

    Parameters
    ----------
    A1_x, A1_y, A1_z : [float]
        Positions.

    A1_m : [float]
        Masses.

    Returns
    -------
    A2_I (3x3 [float])
        The inertia tensor (See Bett et al. (2007), Eqn. (6) etc).
    """
    # Inertia tensor
    A2_I = np.empty((3, 3))

    # Calculate the inertia tensor elements
    for i in range(3):
        for j in range(3):
            A2_I[i, j] = sum(
                A1_m
                * (
                    (A1_x**2 + A1_y**2 + A1_z**2) * delta_ij(i, j)
                    - [A1_x, A1_y, A1_z][i] * [A1_x, A1_y, A1_z][j]
                )
            )

    return A2_I


def Roche_limit(R_planet, rho_planet, rho_moon):
    """Calculate the Roche limit (fluid).

    Parameters
    ----------
    R_planet : float
        Planet radius (kg).

    rho_planet : float
        Planet density (kg m^-3).

    rho_moon : float
        Satellite density (kg m^-3).

    Returns
    -------
    R_Roche : float
        Roche limit (m).
    """
    return 2.455 * R_planet * (rho_planet / rho_moon) ** (1 / 3)


def centre_of_mass(A1_m, A2_pos):
    """Return the centre of mass of a collection of points.

    Parameters
    ----------
    A1_m : [float]
        Masses.

    A2_pos : [[float]]
        Positions ([[x0, y0, z0], [x1, y1, z1], ...]).

    Returns
    -------
    A1_pos_com : [float]
        Centre-of-mass coordinates.
    """
    return sum((np.array(A1_m) * np.array(A2_pos).T).T, axis=0) / sum(A1_m)


def v_esc_from_M_R(M, R):
    """Return the escape velocity from a mass M at a distance R."""
    return np.sqrt(2 * G * M / R)


# ========
# Type checks etc
# ========
def check_bool(param):
    """Check that an input variable is a valid boolean-like string, return the
    appropriate boolean value or raise an error.

    Parameters
    ----------
    param : str
        User-input argument that should be "yes", "no", or similar.

    Returns
    -------
    true_false : bool
        The input converted into the revelant bool.
    """
    # Valid versions of True and False
    A1_true = [True, 1, "True", "true", "1", "yes", "y", "Yes", "Y"]
    A1_false = [False, 0, "False", "false", "0", "no", "n", "No", "N"]

    if param in A1_true:
        true_false = True
    elif param in A1_false:
        true_false = False
    else:
        raise ValueError("Invalid bool input!")

    return true_false


def check_bool_or_def(param, default):
    """Like check_bool() but returns the default if "." or None is the input."""
    if param == "." or param is None:
        return default
    else:
        return check_bool(param)


def check_int(param):
    """Check that an input variable is a valid integer-like string, Return the
    appropriate integer, rounded if necessary, or raise an error.

    Parameters
    ----------
    param : str
        User-input argument that should be an integer as a string.

    Returns
    -------
    param : int
        The input converted into the revelant integer.
    """
    try:
        return int(np.round(float(param)))
    except ValueError:
        raise ValueError("Invalid integer input!")


def check_int_or_def(param, default):
    """Like check_int() but returns the default if "." or None is the input."""
    if param == "." or param is None:
        return default
    else:
        return check_int(param)


def check_float(param):
    """Check that an input variable is a valid float-like string, Return the
    appropriate float, rounded if necessary, or raise an error.

    Parameters
    ----------
    param : str
        User-input argument that should be an float as a string.

    Returns
    -------
    param : float
        The input converted floato the revelant float.
    """
    try:
        return float(param)
    except ValueError:
        raise ValueError("Invalid float input!")


def check_float_or_def(param, default):
    """Like check_float() but returns the default if "." or None is the input."""
    if param == "." or param is None:
        return default
    else:
        return check_float(param)


def check_none(param):
    """Check if an input variable is a valid None-like string, return None if
    so.

    Parameters
    ----------
    param : ?
        User-input argument.

    Returns
    -------
    param : ?
        The input either unchanged or set to None.
    """
    if param in ["None", "none"]:
        return None
    else:
        return param


def check_none_or_def(param, default):
    """Like check_none() but returns the default if "." or None is the input."""
    if param == "." or param is None:
        return default
    else:
        return check_none(param)


def check_option(param, A1_option):
    """Check that a parameter is one of a list of valid options.

    Parameters
    ----------
    param : str
        The parameter to check.

    A1_option : [str]
        A list of valid options.
    """
    if param in A1_option:
        return param
    else:
        message = "\n  Valid options:"
        for option in A1_option:
            message += "\n    %s" % option

        raise ValueError("Invalid input!" + message)


def check_option_or_def(param, A1_option, default):
    """Like check_option() but returns the default if "." or None is the input."""
    if param == "." or param is None:
        return default
    else:
        return check_option(param, A1_option)


def check_end(string, end):
    """Check that a string ends with the required characters and append them
    if not.
    """
    if string[-len(end) :] != end:
        string += end

    return string


def is_between(x, x_min, x_max):
    """Return True if x_min <= x < x_max, False otherwise."""
    if x_min <= x and x < x_max:
        return True
    else:
        return False


def is_within(x, x_mid, dx):
    """Return True if x_mid - dx <= x < x_mid + dx, False otherwise."""
    return is_between(x, x_mid - dx, x_mid + dx)


def where_between(A1_x, x_min, x_max):
    """Return the index array of where x_min <= A1_x < x_max."""
    return np.where((x_min < A1_x) & (A1_x < x_max))[0]


def where_within(A1_x, x_mid, dx):
    """Return the index array of where x_mid - dx <= A1_x < x_mid + dx."""
    return where_between(A1_x, x_mid - dx, x_mid + dx)


def array_from_string(A1_x, type=str):
    """Convert a string user-input array to an actual array, assuming the form
    is e.g. "[x_1,x_2,x_3]" including the brackets and without any spaces.

    Parameters
    ----------
    A1_x : str
        The string form of the array.

    type : type (opt.)
        The desired type of the elements in the converted array.

    Returns
    -------
    A1_x : [type]
        The array.
    """
    assert (
        A1_x[0] == "[" and A1_x[-1] == "]"
    ), "The array must be a string of the form: [x_1,x_2,...]"

    return np.array([type(x) for x in A1_x[1:-1].split(",")])


def ensure_list(A1_x):
    """Convert to a list if not already a list or a numpy array.

    Useful for example in a function that operates over multiple entries in an
    array but also wants to accept a single entry being passed as an argument.
    """
    if not isinstance(A1_x, (list, np.ndarray)):
        return [A1_x]
    else:
        return A1_x


# ========
# Printing utilities
# ========
def auto_format_nice(value):
    """Return a print-ready string of a number with probably-nice formatting."""
    if type(value) is int:
        return "%d" % value
    elif abs(value) < 1e-2:
        return "%.4e" % value
    elif abs(value) < 1e2:
        return "%.5f" % value
    elif abs(value) < 1e4:
        return "%.2f" % value
    else:
        return "%.4e" % value

    def print_mmm(A1_value, middle="mean", prefix="\n", do_header=False):
        """Print the minimum, mean, and maximum of an array for easy checking.

        Parameters
        ----------
        A1_value : [float] or [int]
            The array of values.

        middle : str (opt.)
            The middle value to print. One of "mean", "median", or "both".

        prefix : str (opt.)
            A string to print before the values.

        do_header : bool (opt.)
            If True, then print a header for the values.
        """
        check_option(middle, ["mean", "median", "both"])

        minimum = np.amin(A1_value)
        maximum = np.amax(A1_value)
        mean = np.mean(A1_value, dtype=np.float64)
        median = np.median(A1_value)

        print(prefix, end="")

        space = 13

        if do_header:
            if middle == "mean":
                print(
                    "%s%s%s"
                    % (
                        add_whitespace("Minimum", space),
                        add_whitespace("Mean", space),
                        "Maximum",
                    )
                )
            elif middle == "median":
                print(
                    "%s%s%s"
                    % (
                        add_whitespace("Minimum", space),
                        add_whitespace("Median", space),
                        "Maximum",
                    )
                )
            else:
                print(
                    "%s%s%s%s"
                    % (
                        add_whitespace("Minimum", space),
                        add_whitespace("Mean", space),
                        add_whitespace("Median", space),
                        "Maximum",
                    )
                )

        # Print in a nice format
        if middle == "mean":
            print(
                "%s%s%s"
                % (
                    add_whitespace(auto_format_nice(minimum), space),
                    add_whitespace(auto_format_nice(mean), space),
                    add_whitespace(auto_format_nice(maximum), space),
                )
            )
        elif middle == "median":
            print(
                "%s%s%s"
                % (
                    add_whitespace(auto_format_nice(minimum), space),
                    add_whitespace(auto_format_nice(median), space),
                    add_whitespace(auto_format_nice(maximum), space),
                )
            )
        else:
            print(
                "%s%s%s%s"
                % (
                    add_whitespace(auto_format_nice(minimum), space),
                    add_whitespace(auto_format_nice(mean), space),
                    add_whitespace(auto_format_nice(median), space),
                    add_whitespace(auto_format_nice(maximum), space),
                )
            )


def add_whitespace(string, space):
    """Return a string for aligned printing with adjusted spaces to account for
    the length of the input.

    Example
    -------
    >>> asdf = 123
    >>> qwerty = 456
    >>> print("%s = %d \n""%s = %d" % (
            add_whitespace("asdf", 12), asdf, add_whitespace("qwerty", 12), qwerty
        ))
    asdf         = 123
    qwerty       = 456
    """
    return "%s" % string + " " * (space - len("%s" % string))


def print_dict(
    Di_key_value, key_label=None, value_label=None, sort="key", value_attr=None
):
    """Neatly print a dictionary's contents.

    Parameters
    ----------
    Di_key_value (dictionary)
        The dictionary to print, with string (or string-able) keys.

    key_label, value_label : str (opt.)
        Headers for the key and label columns to print at the start and
        end.

    sort : str (opt.)
        "key":      (Default) Print the items in key order.
        "value":    Print the items in value order.

    value_attr : str (opt.)
        If provided, then instead of printing the values, print a property,
        e.g. for name attributes of custom-class objects.
    """
    check_option(sort, ["key", "value"])

    # Whitespace
    space = max(len("%s" % key) for key in Di_key_value)

    # Header
    if key_label != None and value_label != None:
        space = max(space, len(key_label) + 2)
        print("# %s # %s" % (add_whitespace(key_label, space), value_label))

    # Keys and values
    if sort == "key":
        for key, value in sorted(Di_key_value.items()):
            print("%s : " % (add_whitespace(key, space)), end="")
            if value_attr is None:
                print(value)
            else:
                print(getattr(value, value_attr))
    else:
        for key, value in sorted(Di_key_value.items(), key=lambda k, v: (v, k)):
            print("%s : " % (add_whitespace(key, space)), end="")
            if value_attr is None:
                print(value)
            else:
                print(getattr(value, value_attr))

    # Footer
    if len(Di_key_value) > 20 and key_label != None and value_label != None:
        print("# %s %s" % (add_whitespace(key_label, space), value_label))


def no_latex_string(string):
    """Make a string safe to print, escaping latex special characters.

    Useful for e.g. using a filename with underscores as a figure title.
    """
    # Escape
    for char in ["_", "^"]:
        string = string.replace(char, "\\" + char)

    # Remove
    for char in ["$", "{", "}", "\\left", "\\right", "\\rm"]:
        string = string.replace(char, "")

    # Replace
    string = string.replace("~", " ")

    return string


def s_if_plural(n):
    """Return "s" for pluralising if n > 1, "" otherwise."""
    if n > 1:
        return "s"
    else:
        return ""


def format_array_string(array, format):
    """Return a print-ready string of an array's contents in a given format.

    Parameters
    ----------
        array ([])
            An array of values that can be printed with the given format.

        format : str
            A printing format, e.g. "%d", "%.5g".

            Custom options:
                "string":   Include quotation marks around each string.
                "dorf":     Int or float if decimal places are needed.

    Returns
    -------
        string : str
            The formatted string.
    """
    string = ""

    # Append each element
    # 1D
    if len(np.shape(array)) == 1:
        for x in array:
            if x is None:
                string += "None, "
            else:
                # Custom formats
                if format == "string":
                    string += '"%s", ' % x
                elif format == "dorf":
                    string += "{0}, ".format(str(round(x, 1) if x % 1 else int(x)))
                # Standard formats
                else:
                    string += "%s, " % (format % x)
    # Recursive for higher dimensions
    else:
        for arr in array:
            string += "%s, " % format_array_string(arr, format)

    # Add brackets and remove the trailing comma
    return "[%s]" % string[:-2]


def format_std_form(value, dp=2):
    """Format a string for printing in nice (latex) standard form notation.

    Parameters
    ----------
        value : float
            The value to format.

        dp : int (opt.)
            The number of decimal places to use.
    """
    format = "%%.%df" % dp

    expon = int(np.floor(np.log10(value)))

    coeff = value / 10**expon

    return r"%s$\times 10^{%d}$" % (format % coeff, expon)


def replace_line_in_file(Fp_file, find, new_line):
    """Find the line containing a string in a file and replace the entire line.

    Parameters
    ----------
    Fp_file : str
        The file to edit.

    find : str
        The string to find for identifying the line to replace.

    new_line : str
        The new line text.
    """
    # Read
    with open(Fp_file, "r") as f:
        data = f.readlines()

    # Find and replace
    for i_line, line in enumerate(data):
        if find in line:
            data[i_line] = new_line

    # Write
    with open(Fp_file, "w") as f:
        f.writelines(data)


# ========
# Misc
# ========
def copy_object(self, copy, A1_exclude=[], A1_set_eq=[]):
    """Copy the attributes from one class object to another, if they have not
    been set in the target object already.

    Parameters
    ----------
    self : obj
        The object to potentially copy to.

    copy : obj
        The object to potentially copy from.

    A1_exclude : [str] (opt.)
        An array of strings of the names of attributes to not be copied.

    A1_set_eq : [str] (opt.)
        An array of strings of the names of attributes to be set equal to the
        copy's attribute instead of deep-copying. So future modifications will
        affect both objects.
    """
    assert type(copy) in type(self).mro()

    for var in dir(copy):
        # Skip properties
        if isinstance(getattr(type(self), var, None), property):
            continue

        if not var.startswith("__") and var not in A1_exclude:
            # Copy if the target's attribute either is None or doesn't exist
            do_copy = False
            try:
                if getattr(self, var) is None:
                    do_copy = True
            except AttributeError:
                do_copy = True

            if do_copy:
                # Set equal
                if var in A1_set_eq:
                    # self.var = copy.var
                    exec("self." + var + " = copy." + var + "")

                # Deep copy
                else:
                    # self.var = deepcopy(copy.var)
                    exec("self." + var + " = deepcopy(copy." + var + ")")


def copy_dict(Di_set, Di_copy):
    """Copy the entries from one dictionary to another, if not already set.

    Parameters
    ----------
    Di_set : {}
        The dictionary to potentially copy to.

    Di_copy : {}
        The dictionary to potentially copy from.
    """
    if Di_set is None:
        Di_set = Di_copy
    else:
        for var, value in Di_copy.items():
            if not var in Di_set.keys():
                Di_set[var] = value

    return Di_set


def set_selected_none_attributes(self, A1_exclude=[]):
    """Set a class object's attributes to None if they are "None" or "none".

    Parameters
    ----------
    self : obj
        The object

    A1_exclude : [str]
        An array of strings of the names of attributes to be ignored.
    """
    for var in dir(self):
        # Skip properties
        if isinstance(getattr(type(self), var, None), property):
            continue

        if not var.startswith("__") and var not in A1_exclude:
            value = getattr(self, var)
            if type(value) is str and value in ["None", "none"]:
                exec("self." + var + " = None")


def dict_value_or_default(Di_check, key, default=None):
    """Return a value from a dictionary if the key exists, or the default otherwise."""
    if key in self.Di_check.keys():
        return self.Di_check[key]
    else:
        return default


def remove_substring(string, *substring):
    """Remove one or more substrings from a string, if present."""
    for substring_i in substring:
        string = string.replace(substring_i, "")

    return string


import cProfile, pstats, io


def profile(fnc):
    """A decorator that uses cProfile to profile a function."""

    def inner(*args, **kwargs):

        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()

        s = io.StringIO()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

        return retval

    return inner


# ========
# Planets etc
# ========
class Planet:
    """Class to store properties of planets and similar space things.

    Data from https://ssd.jpl.nasa.gov.

    Parameters
    ----------
    name : str
        The body's name.

    M : float
        Mass (kg).

    R : float
        Mean radius (m).

    R_eq : float (opt.)
        Equatorial radius (m).

    day : float (opt.)
        Day rotation period (s).

    primary : Planet (opt.)
        The primary body around which this orbits.

    a : float (opt.)
        The semi-major axis of the body's orbit (m).

    e : float (opt.)
        The eccentricity of the body's orbit.

    i : float (opt.)
        The inclination of the body's orbit (rad).

    obliquity : float (opt.)
        Axial tilt with respect to the body's orbital plane (rad).

    Attributes
    ----------
    o : orbits.Orbit
        The body's orbit around its primary, if relevant.

    V : float
        Volume (m^-3)

    rho : float
        Density (kg m^-3)

    v_esc : float
        Escape speed (m s^-1)
    """

    def __init__(
        self,
        name,
        M,
        R,
        R_eq=None,
        day=None,
        primary=None,
        a=None,
        e=None,
        i=None,
        obliquity=None,
    ):
        self.name = name
        self.M = M
        self.R = R
        self.R_eq = R_eq
        self.day = day
        self.primary = primary
        self.a = a
        self.e = e
        self.i = i
        self.obliquity = obliquity

        # Derived properties
        self.V = 4 / 3 * pi * R**3
        self.rho = self.M / self.V
        self.v_esc = v_esc_from_M_R(self.M, self.R)

    @property
    def o(self):
        if self.primary is not None:
            return Orbit(
                M_p=self.primary.M,
                m=self.M,
                a=self.a,
                e=self.e,
                i=self.i,
            )
        else:
            return None


if True:
    Su = Planet("Sun", M=1.988475e30, R=6.957e8, day=27 * day_to_s)
    Me = Planet(
        name="Mercury",
        M=0.330103e24,
        R=2439.4e3,
        R_eq=2440.5e3,
        day=58.6463 * day_to_s,
        primary=Su,
        a=0.38709927 * au_to_m,
        e=0.20563593,
        i=7.00497902 * deg_to_rad,
        obliquity=2.1 * arcmin_to_rad,
    )
    Ve = Planet(
        name="Venus",
        M=4.86731e24,
        R=6051.8e3,
        R_eq=6051.8e3,
        day=243.018 * day_to_s,
        primary=Su,
        a=0.72333566 * au_to_m,
        e=0.00677672,
        i=3.39467605 * deg_to_rad,
        obliquity=177.3 * deg_to_rad,
    )
    Ea = Planet(
        name="Earth",
        M=5.9722e24,
        R=6.3710e6,
        R_eq=6378.1366e3,
        day=0.99726968 * day_to_s,
        primary=Su,
        a=1.00000261 * au_to_m,
        e=0.01671123,
        i=-0.00001531 * deg_to_rad,
        obliquity=23.4392911 * deg_to_rad,
    )
    Ma = Planet(
        name="Mars",
        M=0.6417e24,
        R=3389.5e3,
        R_eq=3396.19e3,
        day=1.02595676 * day_to_s,
        primary=Su,
        a=1.52371034 * au_to_m,
        e=0.09339410,
        i=1.84969142 * deg_to_rad,
        obliquity=25.19 * deg_to_rad,
    )
    Ju = Planet(
        name="Jupiter",
        M=1898.125e24,
        R=69911e3,
        R_eq=71492e3,
        day=0.41354 * day_to_s,
        primary=Su,
        a=5.20288700 * au_to_m,
        e=0.04838624,
        i=1.30439695 * deg_to_rad,
        obliquity=3.13 * deg_to_rad,
    )
    Sa = Planet(
        name="Saturn",
        M=568.317e24,
        R=58232e3,
        R_eq=60268e3,
        day=0.44401 * day_to_s,
        primary=Su,
        a=9.53667594 * au_to_m,
        e=0.05386179,
        i=2.48599187 * deg_to_rad,
        obliquity=26.73 * deg_to_rad,
    )
    Ur = Planet(
        name="Uranus",
        M=86.8099e24,
        R=25362e3,
        R_eq=25559e3,
        day=0.71833 * day_to_s,
        primary=Su,
        a=19.18916464 * au_to_m,
        e=0.04725744,
        i=0.77263783 * deg_to_rad,
        obliquity=97.77 * deg_to_rad,
    )
    Ne = Planet(
        name="Neptune",
        M=102.4092e24,
        R=24622e3,
        R_eq=24764e3,
        day=0.67125 * day_to_s,
        primary=Su,
        a=30.06992276 * au_to_m,
        e=0.00859048,
        i=1.77004347 * deg_to_rad,
        obliquity=28.32 * deg_to_rad,
    )

    Mo = Planet(
        name="Moon",
        M=7.342e22,
        R=1.738e6,
        day=6.3872 * day_to_s,
        primary=Ea,
        a=3.844e8,
        e=0.0549,
        i=5.145 * deg_to_rad,
        obliquity=6.67 * deg_to_rad,
    )

    Ph = Planet(
        name="Phobos",
        M=1.07e16,
        R=11.1e3,
        primary=Ma,
        day=0.319 * day_to_s,
        a=9.3772e3,
        e=0.0151,
        i=1.082 * deg_to_rad,
    )
    De = Planet(
        name="Deimos",
        M=1.48e15,
        R=6.3e3,
        primary=Ma,
        day=1.263 * day_to_s,
        a=23.4632e3,
        e=0.00033,
        i=1.791 * deg_to_rad,
    )

    Mi = Planet(
        name="Mimas",
        M=3.75054e19,
        R=197.501e3,
        primary=Ma,
        day=0.942 * day_to_s,
        a=29.17 * Ea.R,
        e=0.0196,
        i=1.572 * deg_to_rad,
    )
    En = Planet(
        name="Enceladus",
        M=1.08097e20,
        R=254.840e3,
        primary=Ma,
        day=1.370 * day_to_s,
        a=37.35 * Ea.R,
        e=0.0047,
        i=0.009 * deg_to_rad,
    )
    Te = Planet(
        name="Tethys",
        M=6.1745e20,
        R=531.1e3,
        primary=Ma,
        day=1.887 * day_to_s,
        a=294619e3,
        e=0.0001,
        i=1.091 * deg_to_rad,
    )
    Di = Planet(
        name="Dione",
        M=1.09530e21,
        R=561.285e3,
        primary=Ma,
        day=2.737 * day_to_s,
        a=56.24 * Ea.R,
        e=0.0022,
        i=0.028 * deg_to_rad,
    )
    Rh = Planet(
        name="Rhea",
        M=2.30646e21,
        R=764.520e3,
        primary=Ma,
        day=4.518 * day_to_s,
        a=82.74 * Ea.R,
        e=0.00126,
        i=0.331 * deg_to_rad,
    )
    Ti = Planet(
        name="Titan",
        M=1.3452e23,
        R=2574.73e3,
        primary=Ma,
        day=15.945 * day_to_s,
        a=1221870e3,
        e=0.0288,
        i=0.28 * deg_to_rad,
    )
    Ia = Planet(
        name="Iapetus",
        M=1.80566e21,
        R=734.4e3,
        primary=Ma,
        day=79.3215 * day_to_s,
        a=3560820e3,
        e=0.02768,
        i=7.489 * deg_to_rad,
    )

    Ves = Planet(
        name="Vesta",
        M=2.589e20,
        R=262.7e3,
        primary=Su,
        day=5.3421 * hour_to_s,
        a=2.36192208 * au_to_m,
        e=0.08944909,
        i=7.14217697 * deg_to_rad,
    )
    Ce = Planet(
        name="Ceres",
        M=938.416e24,
        R=469.7e3,
        R_eq=482.1e3,
        primary=Su,
        day=0.37809 * day_to_s,
        a=2.76725436 * au_to_m,
        e=0.07891253,
        i=10.5868796 * deg_to_rad,
    )
    Pl = Planet(
        name="Pluto",
        M=13029.0e24,
        R=1188.3e3,
        R_eq=1188.3e3,
        primary=Su,
        day=6.3872 * day_to_s,
        a=39.588629 * au_to_m,
        e=0.251837878,
        i=17.14771 * deg_to_rad,
    )

    # Earth units
    Ea_to_SI = Conversions(m=Ea.M, l=Ea.R, t=1)
    SI_to_Ea = Ea_to_SI.inv()

    # Sun/au/year G=1 units
    G1_to_SI = Conversions(m=Su.M, l=au_to_m, t=yr_to_s / (2 * pi))
    SI_to_G1 = G1_to_SI.inv()


# ========
# Colours and lines
# ========
if True:
    # List of nice, contrasting colours for general use
    # The first few are well-separated in luminosity for greyscale too
    A1_c = [
        "#1199ff",
        "#ee4400",
        "#7711dd",
        "#44dd44",
        "#ffdd00",
        "#60fff0",
        "#775533",
        "#ff77dd",
        "#707070",
        "#ccaaee",
        "#0000cc",
        "#cc0000",
        "#660077",
        "#007700",
        "#ff9922",
        "#aaccff",
        "#55ffbb",
        "#ffcccc",
        "#cc8866",
        "#111111",
        "#aaaaaa",
        "#660000",
        "#9922cc",
        "#004411",
    ]
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=A1_c)

    # ========
    # Tweaked colour maps
    # ========
    # Remove the dark low end of inferno
    cmap = plt.get_cmap("inferno")
    new_colours = cmap(np.linspace(0.15, 1, cmap.N))
    cmap_melty = mpl.colors.LinearSegmentedColormap.from_list("melty", new_colours)
    cmap_melty_r = mpl.colors.LinearSegmentedColormap.from_list(
        "melty_r", cmap_melty(np.linspace(1, 0, cmap_melty.N))
    )

    # Remove the dark low end of viridis
    cmap = plt.get_cmap("viridis")
    new_colours = cmap(np.linspace(0.15, 1, cmap.N))
    cmap_viridis_trim = mpl.colors.LinearSegmentedColormap.from_list(
        "viridis_trim", new_colours
    )
    cmap_viridis_trim_r = mpl.colors.LinearSegmentedColormap.from_list(
        "viridis_trim_r", cmap_viridis_trim(np.linspace(1, 0, cmap_viridis_trim.N))
    )

    # Black opacity
    new_colours = np.zeros((256, 4))
    new_colours[:, 3] = 1 - np.tanh(np.arange(256) / 150)
    cmap_k_opac_tanh = mpl.colors.LinearSegmentedColormap.from_list(
        "k_opac_tanh", new_colours
    )
    new_colours[:, 3] = np.tanh(np.arange(256) / 150)
    cmap_k_opac_tanh_r = mpl.colors.LinearSegmentedColormap.from_list(
        "k_opac_tanh_r", new_colours
    )
    new_colours[:, 3] = np.linspace(1, 0, 256)
    cmap_k_opac_lin = mpl.colors.LinearSegmentedColormap.from_list(
        "k_opac_lin", new_colours
    )
    new_colours[:, 3] = np.linspace(0, 1, 256)
    cmap_k_opac_lin_r = mpl.colors.LinearSegmentedColormap.from_list(
        "k_opac_lin_r", new_colours
    )

    # White opacity
    new_colours = np.ones((256, 4))
    new_colours[:, 3] = 1 - np.tanh(np.arange(256) / 150)
    cmap_w_opac_tanh = mpl.colors.LinearSegmentedColormap.from_list(
        "w_opac_tanh", new_colours
    )
    new_colours[:, 3] = np.tanh(np.arange(256) / 150)
    cmap_w_opac_tanh_r = mpl.colors.LinearSegmentedColormap.from_list(
        "w_opac_tanh_r", new_colours
    )
    new_colours[:, 3] = np.linspace(1, 0, 256)
    cmap_w_opac_lin = mpl.colors.LinearSegmentedColormap.from_list(
        "w_opac_lin", new_colours
    )
    new_colours[:, 3] = np.linspace(0, 1, 256)
    cmap_w_opac_lin_r = mpl.colors.LinearSegmentedColormap.from_list(
        "w_opac_lin_r", new_colours
    )

    # Other perceptually uniform colour maps
    path_cmap = path + "/tables"
    os.makedirs(path_cmap, exist_ok=True)
    try:
        # Get the maps in an environment with colorcet installed and save them for
        # later use elsewhere
        import colorcet

        cmap_rbow = colorcet.m_rainbow
        cmap_rbow2 = colorcet.m_rainbow_bgyrm_35_85_c69
        cmap_rbowj = colorcet.m_rainbow_bgyr_10_90_c83
        cmap_bwr = colorcet.m_diverging_bwr_40_95_c42
        cmap_gwv = colorcet.m_diverging_gwv_55_95_c39
        cmap_cwr = colorcet.m_diverging_tritanopic_cwr_75_98_c20
        cmap_cgo = colorcet.m_isoluminant_cgo_80_c38
        cmap_cjo = colorcet.m_diverging_isoluminant_cjo_70_c25
        cmap_cm = colorcet.m_isoluminant_cm_70_c39
        with open(path_cmap + "/cmap_rbow.pkl", "wb") as f:
            pickle.dump(cmap_rbow, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_rbow2.pkl", "wb") as f:
            pickle.dump(cmap_rbow2, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_rbowj.pkl", "wb") as f:
            pickle.dump(cmap_rbowj, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_bwr.pkl", "wb") as f:
            pickle.dump(cmap_bwr, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_gwv.pkl", "wb") as f:
            pickle.dump(cmap_gwv, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_cwr.pkl", "wb") as f:
            pickle.dump(cmap_cwr, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_cgo.pkl", "wb") as f:
            pickle.dump(cmap_cgo, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_cjo.pkl", "wb") as f:
            pickle.dump(cmap_cjo, f, pickle.HIGHEST_PROTOCOL)
        with open(path_cmap + "/cmap_cm.pkl", "wb") as f:
            pickle.dump(cmap_cm, f, pickle.HIGHEST_PROTOCOL)
    except ImportError:
        # Load the previously stored maps
        with open(path_cmap + "/cmap_rbow.pkl", "rb") as f:
            cmap_rbow = pickle.load(f)
        with open(path_cmap + "/cmap_rbow2.pkl", "rb") as f:
            cmap_rbow2 = pickle.load(f)
        with open(path_cmap + "/cmap_rbowj.pkl", "rb") as f:
            cmap_rbowj = pickle.load(f)
        with open(path_cmap + "/cmap_bwr.pkl", "rb") as f:
            cmap_bwr = pickle.load(f)
        with open(path_cmap + "/cmap_gwv.pkl", "rb") as f:
            cmap_gwv = pickle.load(f)
        with open(path_cmap + "/cmap_cwr.pkl", "rb") as f:
            cmap_cwr = pickle.load(f)
        with open(path_cmap + "/cmap_cgo.pkl", "rb") as f:
            cmap_cgo = pickle.load(f)
        with open(path_cmap + "/cmap_cjo.pkl", "rb") as f:
            cmap_cjo = pickle.load(f)
        with open(path_cmap + "/cmap_cm.pkl", "rb") as f:
            cmap_cm = pickle.load(f)
    cmap_rbow_r = mpl.colors.LinearSegmentedColormap.from_list(
        "rbow_r", cmap_rbow(np.linspace(1, 0, cmap_rbow.N))
    )
    cmap_rbow2_r = mpl.colors.LinearSegmentedColormap.from_list(
        "rbow2_r", cmap_rbow2(np.linspace(1, 0, cmap_rbow2.N))
    )
    cmap_rbowj_r = mpl.colors.LinearSegmentedColormap.from_list(
        "rbowj_r", cmap_rbowj(np.linspace(1, 0, cmap_rbowj.N))
    )
    cmap_bwr_r = mpl.colors.LinearSegmentedColormap.from_list(
        "bwr_r", cmap_bwr(np.linspace(1, 0, cmap_bwr.N))
    )
    cmap_gwv_r = mpl.colors.LinearSegmentedColormap.from_list(
        "gwv_r", cmap_gwv(np.linspace(1, 0, cmap_gwv.N))
    )
    cmap_cwr_r = mpl.colors.LinearSegmentedColormap.from_list(
        "cwr_r", cmap_cwr(np.linspace(1, 0, cmap_cwr.N))
    )
    cmap_cgo_r = mpl.colors.LinearSegmentedColormap.from_list(
        "cgo_r", cmap_cgo(np.linspace(1, 0, cmap_cgo.N))
    )
    cmap_cjo_r = mpl.colors.LinearSegmentedColormap.from_list(
        "cjo_r", cmap_cjo(np.linspace(1, 0, cmap_cjo.N))
    )
    cmap_cm_r = mpl.colors.LinearSegmentedColormap.from_list(
        "cm_r", cmap_cm(np.linspace(1, 0, cmap_cm.N))
    )

    # Matplotlib named colours and hex codes
    Di_c_hex = {
        "aliceblue": "#F0F8FF",
        "antiquewhite": "#FAEBD7",
        "aqua": "#00FFFF",
        "aquamarine": "#7FFFD4",
        "azure": "#F0FFFF",
        "beige": "#F5F5DC",
        "bisque": "#FFE4C4",
        "black": "#000000",
        "blanchedalmond": "#FFEBCD",
        "blue": "#0000FF",
        "blueviolet": "#8A2BE2",
        "brown": "#A52A2A",
        "burlywood": "#DEB887",
        "cadetblue": "#5F9EA0",
        "chartreuse": "#7FFF00",
        "chocolate": "#D2691E",
        "coral": "#FF7F50",
        "cornflowerblue": "#6495ED",
        "cornsilk": "#FFF8DC",
        "crimson": "#DC143C",
        "cyan": "#00FFFF",
        "darkblue": "#00008B",
        "darkcyan": "#008B8B",
        "darkgoldenrod": "#B8860B",
        "darkgray": "#A9A9A9",
        "darkgreen": "#006400",
        "darkkhaki": "#BDB76B",
        "darkmagenta": "#8B008B",
        "darkolivegreen": "#556B2F",
        "darkorange": "#FF8C00",
        "darkorchid": "#9932CC",
        "darkred": "#8B0000",
        "darksalmon": "#E9967A",
        "darkseagreen": "#8FBC8F",
        "darkslateblue": "#483D8B",
        "darkslategray": "#2F4F4F",
        "darkturquoise": "#00CED1",
        "darkviolet": "#9400D3",
        "deeppink": "#FF1493",
        "deepskyblue": "#00BFFF",
        "dimgray": "#696969",
        "dodgerblue": "#1E90FF",
        "firebrick": "#B22222",
        "floralwhite": "#FFFAF0",
        "forestgreen": "#228B22",
        "fuchsia": "#FF00FF",
        "gainsboro": "#DCDCDC",
        "ghostwhite": "#F8F8FF",
        "gold": "#FFD700",
        "goldenrod": "#DAA520",
        "gray": "#808080",
        "green": "#008000",
        "greenyellow": "#ADFF2F",
        "honeydew": "#F0FFF0",
        "hotpink": "#FF69B4",
        "indianred": "#CD5C5C",
        "indigo": "#4B0082",
        "ivory": "#FFFFF0",
        "khaki": "#F0E68C",
        "lavender": "#E6E6FA",
        "lavenderblush": "#FFF0F5",
        "lawngreen": "#7CFC00",
        "lemonchiffon": "#FFFACD",
        "lightblue": "#ADD8E6",
        "lightcoral": "#F08080",
        "lightcyan": "#E0FFFF",
        "lightgoldenrodyellow": "#FAFAD2",
        "lightgreen": "#90EE90",
        "lightgray": "#D3D3D3",
        "lightpink": "#FFB6C1",
        "lightsalmon": "#FFA07A",
        "lightseagreen": "#20B2AA",
        "lightskyblue": "#87CEFA",
        "lightslategray": "#778899",
        "lightsteelblue": "#B0C4DE",
        "lightyellow": "#FFFFE0",
        "lime": "#00FF00",
        "limegreen": "#32CD32",
        "linen": "#FAF0E6",
        "magenta": "#FF00FF",
        "maroon": "#800000",
        "mediumaquamarine": "#66CDAA",
        "mediumblue": "#0000CD",
        "mediumorchid": "#BA55D3",
        "mediumpurple": "#9370DB",
        "mediumseagreen": "#3CB371",
        "mediumslateblue": "#7B68EE",
        "mediumspringgreen": "#00FA9A",
        "mediumturquoise": "#48D1CC",
        "mediumvioletred": "#C71585",
        "midnightblue": "#191970",
        "mintcream": "#F5FFFA",
        "mistyrose": "#FFE4E1",
        "moccasin": "#FFE4B5",
        "navajowhite": "#FFDEAD",
        "navy": "#000080",
        "oldlace": "#FDF5E6",
        "olive": "#808000",
        "olivedrab": "#6B8E23",
        "orange": "#FFA500",
        "orangered": "#FF4500",
        "orchid": "#DA70D6",
        "palegoldenrod": "#EEE8AA",
        "palegreen": "#98FB98",
        "paleturquoise": "#AFEEEE",
        "palevioletred": "#DB7093",
        "papayawhip": "#FFEFD5",
        "peachpuff": "#FFDAB9",
        "peru": "#CD853F",
        "pink": "#FFC0CB",
        "plum": "#DDA0DD",
        "powderblue": "#B0E0E6",
        "purple": "#800080",
        "red": "#FF0000",
        "rosybrown": "#BC8F8F",
        "royalblue": "#4169E1",
        "saddlebrown": "#8B4513",
        "salmon": "#FA8072",
        "sandybrown": "#FAA460",
        "seagreen": "#2E8B57",
        "seashell": "#FFF5EE",
        "sienna": "#A0522D",
        "silver": "#C0C0C0",
        "skyblue": "#87CEEB",
        "slateblue": "#6A5ACD",
        "slategray": "#708090",
        "snow": "#FFFAFA",
        "springgreen": "#00FF7F",
        "steelblue": "#4682B4",
        "tan": "#D2B48C",
        "teal": "#008080",
        "thistle": "#D8BFD8",
        "tomato": "#FF6347",
        "turquoise": "#40E0D0",
        "violet": "#EE82EE",
        "wheat": "#F5DEB3",
        "white": "#FFFFFF",
        "whitesmoke": "#F5F5F5",
        "yellow": "#FFFF00",
        "yellowgreen": "#9ACD32",
    }

    # ========
    # Line styles
    # ========
    ls_dot = (0, (1, 4))
    ls_dot_dense = (0, (1, 2))
    ls_dot_space = (0, (1, 6))
    ls_dash = (0, (6, 5))
    ls_dash_dot = (0, (6, 3, 1, 3))
    ls_dash_long = (0, (8, 8))
    ls_dash_longer = (0, (16, 16))
    ls_dash_dot_dot = (0, (7, 3, 1, 3, 1, 3))
    ls_dash_dash_dot = (0, (6, 1, 6, 2, 1, 2))
    ls_dash_dash_dot_dot = (0, (5, 1, 5, 2, 1, 1, 1, 2))
    ls_dash_dot_dot_dot = (0, (5, 2, 1, 1, 1, 1, 1, 2))
    ls_dash_dash_dash_dot = (0, (4, 1, 4, 1, 4, 2, 1, 2))
    ls_dash_dot_dot_dot_dot = (0, (8, 1, 1, 1, 1, 1, 1, 1, 1, 1))

    A1_ls = [
        "-",
        ls_dash,
        ls_dot,
        ls_dash_dot,
        ls_dash_dot_dot,
        ls_dash_dash_dot,
        ls_dash_dash_dot_dot,
        ls_dash_dot_dot_dot,
        ls_dash_dash_dash_dot,
        ls_dash_dot_dot_dot_dot,
    ]


class MidNorm(mpl.colors.Normalize):
    """Normalise a colour bar to centre on a chosen midpoint."""

    def __init__(self, vcenter=None, vmin=None, vmax=None, clip=False):
        self.vcenter = vcenter
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def hex_to_grey(hex):
    """Convert a hexadecimal string to its greyscale luminosity.

    Parameters
    ----------
    hex : str
        Hexadecimal colour: "#rrggbb".

    Returns
    -------
    luminosity : float
        The weighted greyscale luminosity.
    """
    r = int(hex[1:3], 16)
    g = int(hex[3:5], 16)
    b = int(hex[5:], 16)

    return 0.21 * r + 0.72 * g + 0.07 * b


def tweak_luminosity(c, tweak):
    """Make a colour ligher or darker.

    Parameters
    ----------
    c : rgb, rgba tuple, or hex string
        The colour to modify.

    tweak : float
        The amount to change the luminosity [0,1], e.g.

        0.1     Much darker
        0.4     Slightly darker
        0.6     Slightly lighter
        0.9     Much lighter

    Returns
    -------
    c : rgba tuple
        The modified colour.

    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    tweak = (1 - tweak) * 2

    if c == "none":
        return c

    # Convert hex to rgb
    if c[0] == "#":
        c = tuple(int(c[i : i + 2], 16) / 255 for i in (1, 3, 5))
        alpha = 1
    # Transparency
    else:
        try:
            alpha = float(c[3])
        except IndexError:
            alpha = 1

    # Convert to HLS
    c = colorsys.rgb_to_hls(*mpl.colors.to_rgb(c))

    # New luminosity
    l = 1 - tweak * (1 - c[1])
    if l < 0:
        l = 0
    elif 1 < l:
        l = 1

    # Back to RGB
    c = colorsys.hls_to_rgb(c[0], l, c[2])

    # Keep the original alpha
    return (c[0], c[1], c[2], alpha)


def A1_colour_from_cmap(cmap, n):
    """Return n equally spaced colours from a colour map."""
    return [cmap(np.linspace(0, 1, n))[i] for i in range(n)]


# ========
# Matplotlib defaults
# ========
font_size = 26
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
params = {
    "backend": "ps",
    "text.latex.preamble": r"\usepackage{gensymb}" r"\usepackage{wasysym}",
    "axes.labelsize": font_size,
    "axes.titlesize": font_size,
    "font.size": font_size,
    "legend.fontsize": font_size - 4,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "text.usetex": True,
    "figure.figsize": [9, 9],
    "font.family": "serif",
    # "font.family"           : "sans-serif",
    # "font.sans-serif"       : "Open Sans",
    "savefig.dpi": 100,
    "legend.framealpha": 0.8,
    "legend.labelspacing": 0.3,
    "lines.linewidth": 1.7,
}
mpl.rcParams.update(params)
# mpl.rc("mathtext", fontset="stixsans")


# ========
# Axes etc
# ========
def set_std_form_axes(ax, force_x=False, force_y=False):
    """Set the axis labels to use standard form, if the current limits are very
    large or small (and the scale is linear).

    Parameters
    ----------
    ax : Axes
        The plot axes.

    force_x, force_y : bool (opt.)
        False   Only set stardard form if the axis limits are very large or small.
        True    Set standard form regardless of the current limits.
    """
    # Check if either axis is logarithmic
    if ax.get_xscale() == "log":
        log_x = True
    else:
        log_x = False
    if ax.get_yscale() == "log":
        log_y = True
    else:
        log_y = False

    # Boundary labels
    lim_min = 1e-3
    lim_max = 1e4

    # ========
    # Current smallest and largest axis labels
    # ========
    A1_x_tick = abs(ax.get_xticks())
    A1_y_tick = abs(ax.get_yticks())
    # Ignore 0
    A1_x_tick = A1_x_tick[A1_x_tick != 0]
    A1_y_tick = A1_y_tick[A1_y_tick != 0]

    # ========
    # Set to standard form
    # ========
    # x axis
    if len(A1_x_tick) > 0:
        x_min = np.amin(A1_x_tick)
        x_max = np.amax(A1_x_tick)
        if (x_min < lim_min or lim_max < x_max or force_x) and not log_x:
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            ax.xaxis.set_major_formatter(formatter)
            try:
                ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
            except AttributeError:
                pass
    # y axis
    if len(A1_y_tick) > 0:
        y_min = np.amin(A1_y_tick)
        y_max = np.amax(A1_y_tick)
        if (y_min < lim_min or lim_max < y_max or force_y) and not log_y:
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            ax.yaxis.set_major_formatter(formatter)
            try:
                ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
            except AttributeError:
                pass


def set_std_form_cbar(cbar, A1_colour, force=False):
    """Set a colour bar's axis labels to use standard form, if the current
    limits are very large or small.

    Parameters
    ----------
    cbar : Colorbar
        The colour bar object.

    A1_colour : [float]
        Array of floats that are setting the colour values.

    force : bool (opt.)
        False   Only set stardard form if the axis limits are very large or small.
        True    Set standard form regardless of the current limits.
    """
    # Boundary limits
    lim_min = 1e-2
    lim_max = 1e3

    # Current axis limits
    c_min = abs(np.amin(A1_colour))
    c_max = abs(np.amax(A1_colour))

    # Set to standard form
    if (c_min < lim_min and c_min != 0) or lim_max < c_max or force:
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()


def set_large_ticks(ax):
    """Set larger ticks on plot axes, especially for logarithmic scales.

    Parameters
    ----------
    ax : Axes
        The plot axes.
    """
    # Check if either axis is logarithmic
    log_x = False
    log_y = False
    if ax.get_xscale() == "log":
        log_x = True
    if ax.get_yscale() == "log":
        log_y = True

    # Tick sizes
    width_log = 1.0
    width_lin = 0.8
    width_min = 0.8
    length_log = 10
    length_lin = 6
    length_min = 4

    # First set all major and minor ticks to the logarithmic sizes
    ax.tick_params("both", width=width_log, length=length_log, which="major")
    ax.tick_params("both", width=width_log, length=length_lin, which="minor")

    # Reset linear ticks ((this weird order seems to work best))
    if not log_x:
        ax.xaxis.set_tick_params(width=width_lin, length=length_lin, which="major")
        ax.xaxis.set_tick_params(width=width_min, length=length_min, which="minor")
    if not log_y:
        ax.yaxis.set_tick_params(width=width_lin, length=length_lin, which="major")
        ax.yaxis.set_tick_params(width=width_min, length=length_min, which="minor")


def set_auto_limits(
    ax,
    A1_x,
    A1_y,
    frac_x=0.05,
    frac_y=0.05,
    set_x=True,
    set_y=True,
    min_x=None,
    min_y=None,
    max_x=None,
    max_y=None,
):
    """Set the axes limits to a nice small distance away from the highest and
    lowest points.

    If a logarithmic axis covers less than two decades, add standard form
    axis labels for some minor ticks as well.

    Parameters
    ----------
    ax : Axes
        The plot axes.

    A1_x, A1_y : [float]
        The data point x and y values.

    frac_x, frac_y : float
        The fraction of the data's range by which to extend the axes.
        Default 1/20.

    set_x, set_y : bool (opt.)
        True    Set new limits for that axis.
        False   Do not set new limits for that axis.

    min_x, min_y : bool (opt.)
        If not None then force the minimum value for that axis.

    max_x, max_y : float (opt.)
        If not None then force the maximum value for that axis.
    """
    # ========
    # Check args
    # ========
    set_x = check_bool(set_x)
    set_y = check_bool(set_y)
    min_x = check_none(min_x)
    min_y = check_none(min_y)
    max_x = check_none(max_x)
    max_y = check_none(max_y)
    A1_x = np.array(A1_x)
    A1_y = np.array(A1_y)

    # Check if either axis is logarithmic
    log_x = False
    log_y = False
    if ax.get_xscale() == "log":
        log_x = True
    if ax.get_yscale() == "log":
        log_y = True

    # ========
    # Axis limits
    # ========
    x_min = np.nanmin(A1_x)
    x_max = np.nanmax(A1_x)
    y_min = np.nanmin(A1_y)
    y_max = np.nanmax(A1_y)
    if log_x:
        # Use the lowest positive value for the log scale minimum
        if x_min <= 0:
            x_min = np.nanmin(A1_x[A1_x > 0])

        x_min = np.log10(x_min)
        x_max = np.log10(x_max)
    if log_y:
        # Use the lowest positive value for the log scale minimum
        if y_min <= 0:
            y_min = np.nanmin(A1_y[A1_y > 0])

        y_min = np.log10(y_min)
        y_max = np.log10(y_max)
    if min_x is not None:
        x_min = min_x
    if min_y is not None:
        y_min = min_y
    if max_x is not None:
        x_max = max_x
    if max_y is not None:
        y_max = max_y

    # Set the limits to be a fraction of the total range outside the points
    x_1 = x_min - (x_max - x_min) * frac_x
    x_2 = x_max + (x_max - x_min) * frac_x
    y_1 = y_min - (y_max - y_min) * frac_y
    y_2 = y_max + (y_max - y_min) * frac_y
    if log_x:
        x_1 = 10 ** (x_min - (x_max - x_min) * frac_x)
        x_2 = 10 ** (x_max + (x_max - x_min) * frac_x)
    if log_y:
        y_1 = 10 ** (y_min - (y_max - y_min) * frac_y)
        y_2 = 10 ** (y_max + (y_max - y_min) * frac_y)

    if min_x is not None:
        x_1 = min_x
    if min_y is not None:
        y_1 = min_y
    if max_x is not None:
        x_2 = max_x
    if max_y is not None:
        y_2 = max_y

    if set_x:
        ax.set_xlim(x_1, x_2)
    if set_y:
        ax.set_ylim(y_1, y_2)

    # ========
    # Logarithmic labels
    # ========
    if log_x and set_x:
        if np.log10(x_2) - np.log10(x_1) < 2 and np.log10(x_2) - np.log10(x_1) > 0.1:
            # Visible ticks
            dec = 10 ** (np.floor(np.log10(x_1)))
            A1_x_tick = (
                np.append(
                    np.arange(np.floor(x_1 / dec) + 1, 10, 1),
                    np.arange(1, np.floor(x_2 / (dec * 10)), 1) * 10,
                )
            ) * dec

            # Tick first digits and order of magnitude
            A1_x_tick_digit = np.array([int(str(tick)[0]) for tick in A1_x_tick])
            A1_x_tick_dec = np.log10(A1_x_tick / A1_x_tick_digit)

            try:
                # Standard form labels
                A1_x_ticklabel = np.array(
                    [
                        r"%d$\times$10$^%d$"
                        % (A1_x_tick_digit[i_tick], A1_x_tick_dec[i_tick])
                        for i_tick in range(len(A1_x_tick))
                    ]
                )

                # No labels other than 5 or 3 and 5
                if np.log10(x_2) - np.log10(x_1) < 1:
                    sel_ticks = np.where(
                        (A1_x_tick / (10 ** np.floor(np.log10(A1_x_tick))) != 3)
                        & (A1_x_tick / (10 ** np.floor(np.log10(A1_x_tick))) != 5)
                    )[0]
                else:
                    sel_ticks = np.where(
                        A1_x_tick / (10 ** np.floor(np.log10(A1_x_tick))) != 5
                    )[0]

                A1_x_ticklabel[sel_ticks] = ""

                ax.xaxis.set_ticks(A1_x_tick, minor=True)
                ax.xaxis.set_ticklabels(A1_x_ticklabel, minor=True)
            except OverflowError:
                ## Arises from x_tick_digit = 0, but why?
                pass
    if log_y and set_y:
        if np.log10(y_2) - np.log10(y_1) < 2 and np.log10(y_2) - np.log10(y_1) > 0.1:
            # Visible ticks
            dec = 10 ** (np.floor(np.log10(y_1)))
            A1_y_tick = (
                np.append(
                    np.arange(np.floor(y_1 / dec) + 1, 10, 1),
                    np.arange(1, np.floor(y_2 / (dec * 10)), 1) * 10,
                )
            ) * dec

            # Tick first digits and order of magnitude
            A1_y_tick_digit = np.array([int(str(tick)[0]) for tick in A1_y_tick])
            A1_y_tick_dec = np.log10(A1_y_tick / A1_y_tick_digit)

            # Standard form labels
            try:
                A1_y_ticklabel = np.array(
                    [
                        r"%d$\times$10$^%d$"
                        % (A1_y_tick_digit[i_tick], A1_y_tick_dec[i_tick])
                        for i_tick in range(len(A1_y_tick))
                    ]
                )

                if np.log10(y_2) - np.log10(y_1) < 1:
                    sel_ticks = np.where(
                        (A1_y_tick / (10 ** np.floor(np.log10(A1_y_tick))) != 3)
                        & (A1_y_tick / (10 ** np.floor(np.log10(A1_y_tick))) != 5)
                    )[0]
                else:
                    sel_ticks = np.where(
                        A1_y_tick / (10 ** np.floor(np.log10(A1_y_tick))) != 5
                    )[0]

                A1_y_ticklabel[sel_ticks] = ""

                ax.yaxis.set_ticks(A1_y_tick, minor=True)
                ax.yaxis.set_ticklabels(A1_y_ticklabel, minor=True)
            except OverflowError:
                pass


def auto_tick_spacing(ax_size):
    """Find a nice spacing for axis ticks.

    Parameters
    ----------
    ax_size : float
        The size of the axis.

    Returns
    -------
    tick_base, tick_base_minor : float
        The spacing for the major and minor ticks.
    """
    # Space for at least three ticks
    size_max = 0.3 * ax_size

    # Order of magnitude
    log10 = np.floor(np.log10(size_max))
    magn = 10**log10

    # Prefactor of 1, 2, or 5, and 0.5 or 1 for minor ticks
    if size_max / magn < 2:
        tick_base = magn
        tick_base_minor = 0.5 * magn
    elif size_max / magn < 5:
        tick_base = 2 * magn
        tick_base_minor = magn
    else:
        tick_base = 5 * magn
        tick_base_minor = magn

    return tick_base, tick_base_minor


def set_cbar_font_size(cbar, fontsize=24):
    """Set the font size for colour bar text.

    Parameters
    ----------
    cbar : Colorbar
        A colour bar object.

    fontsize : float (opt.)
        The font size.
    """
    cbar.ax.tick_params(labelsize=fontsize - 2)
    cbar.ax.yaxis.label.set_font_properties(
        mpl.font_manager.FontProperties(size=fontsize)
    )


def set_font_size(ax, cbar=None, fontsize=24):
    """Set the font size for all plot text.

    Parameters
    ----------
    ax : Axes
        The plot axes.

    cbar : Colorbar (opt.)
        A colour bar object.

    fontsize : float (opt.)
        The font size.
    """
    # Main labels
    for item in (
        [ax.title, ax.xaxis.label, ax.yaxis.label]
        + ax.get_xticklabels()
        + ax.get_yticklabels()
    ):
        item.set_fontsize(fontsize)

    # Minor tick labels
    for item in [] + ax.xaxis.get_minorticklabels() + ax.yaxis.get_minorticklabels():
        item.set_fontsize(fontsize - 2)

    # Colour bar
    if cbar:
        set_cbar_font_size(cbar, fontsize=fontsize)


def add_far_ticks(ax, do_x=True, do_y=True):
    """Add axis ticks on the top and/or right axes."""
    if do_x:
        ax_t = ax.secondary_xaxis("top")
        ax_t.set_xticks(ax.get_xticks())

        # Log or minor ticks
        if ax.get_xscale() == "log":
            ax_t.set_xscale("log")
        else:
            A1_minor_ticks = ax.get_xticks(minor=True)
            if len(A1_minor_ticks) > 0:
                ax_t.set_xticks(A1_minor_ticks, minor=True)

        # Inward ticks, no labels
        ax_t.tick_params(axis="x", which="both", direction="in")
        ax_t.set_xticklabels([])
        ax_t.set_xticklabels([], minor=True)
        set_large_ticks(ax_t)
    else:
        ax_t = None

    if do_y:
        ax_r = ax.secondary_yaxis("right")
        ax_r.set_yticks(ax.get_yticks())

        # Log or minor ticks
        if ax.get_yscale() == "log":
            ax_r.set_yscale("log")
        else:
            A1_minor_ticks = ax.get_yticks(minor=True)
            if len(A1_minor_ticks) > 0:
                ax_r.set_yticks(A1_minor_ticks, minor=True)

        # Inward ticks, no labels
        ax_r.tick_params(axis="y", which="both", direction="in")
        ax_r.set_yticklabels([])
        ax_r.set_yticklabels([], minor=True)
        set_large_ticks(ax_r)
    else:
        ax_r = None

    return ax_t, ax_r


def tight_layout():
    """Call plt.tight_layout() skipping warnings and avoiding lockfile errors."""
    # Skip warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # # Try to avoid lockfile errors
        # for i in range(10):
        #     try:
        #         plt.tight_layout()
        #     except (TimeoutError, FileNotFoundError):
        #         sleep(1 + 3 * np.random.random())

        plt.tight_layout()


def savefig(*args, **kwargs):
    """Call plt.savefig() avoiding lockfile errors."""
    # Make the directory if it doesn't exist
    os.makedirs(os.path.dirname(args[0]), exist_ok=True)

    # # Try to avoid lockfile errors
    # for i in range(10):
    #     try:
    #         plt.savefig(*args, **kwargs)
    #     except (FileNotFoundError, TimeoutError):
    #         sleep(1 + 3 * np.random.random())

    plt.savefig(*args, **kwargs)


def auto_nice_plot(A1_ax=None, cbar=None, fontsize=24, do_far_ticks=True):
    """Call standard nice-plot-making functions.

    Parameters
    ----------
    A1_ax : Axes or [plt.Axes]
        One or more plot axes.
    """
    if A1_ax is None:
        A1_ax = plt.gcf().axes

    for ax in ensure_list(A1_ax):
        set_large_ticks(ax)
        set_std_form_axes(ax)
        set_font_size(ax, fontsize=fontsize)
        if do_far_ticks:
            add_far_ticks(ax)

    if cbar is not None:
        set_cbar_font_size(cbar, fontsize=fontsize)

    tight_layout()


# ========
# Plotting functions
# ========
def plot_lines(
    A1_x,
    A1_y,
    c="k",
    A1_p_c=None,
    cmap=None,
    vmin=None,
    vmax=None,
    norm=mpl.colors.Normalize,
    lw=1,
    A1_lw=None,
    ls="-",
    ls_step_on=None,
    ls_step_off=None,
    A1_ax_lim=None,
    alpha=1,
    zorder=None,
    label=None,
    ax=None,
):
    """Plot a line collection with optional gradients, e.g. by colour.

    Parameters
    ----------
    A1_x, A1_y : [float]
        The x and y coordinates to plot.

    c, A1_p_c : str, [float] (opt.)
        The colour and/or the colour parameter values for each segment.

    cmap  : Colormap or str (opt.)
        The colour map to use. Or:
            fade                Tweak alpha.
            tweak_luminosity    Tweak the luminosity.

    vmin, vmax, norm : float, float, Normalize (opt.)
        The min and max values and the normalisation for the colour map.

    lw or A1_lw : float or [float] (opt.)
        The linewidth or the linewidths for each segment.

    ls : str
        The linestyle, done as a simple workaround, choose from: "-", "--", ":".
        Requires axis limits to be set or provided to be reliably effective.

    ls_step_on, ls_step_off : float, float (opt.)
        Instead of ls, set the linestyle directly with the lengths of the on and
        off segments, as fractions of the axis side.

    A1_ax_lim : [float] (opt.)
        The axis limits [x_min, x_max, y_min, y_max]. Allows ls options to work
        well, and avoids plotting outside the limits. Set to `"auto"` to use the
        figure values.

    alpha, zorder, label : str, str, float, float, float, str (opt.)
        The opacity, zorder, and label.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.
    """
    if ax is None:
        ax = plt.gca()
    if ls_step_on is not None:
        ls = None

    # Array of single segments [ ..., [ [ x_i, y_i ], [ x_i+1, y_i+1 ] ], ... ]
    num_seg = len(A1_x) - 1
    segments = np.empty((num_seg, 2, 2))
    segments[:, 0, 0] = A1_x[:-1]
    segments[:, 1, 0] = A1_x[1:]
    segments[:, 0, 1] = A1_y[:-1]
    segments[:, 1, 1] = A1_y[1:]

    # Colours
    # Base
    A1_colour = np.array([list(mpl.colors.to_rgba(c))] * num_seg)
    # By parameter
    if A1_p_c is not None:
        # Normalise the colour parameter values to between 0 and 1
        if vmin is None:
            vmin = np.amin(A1_p_c)
        if vmax is None:
            vmax = np.amax(A1_p_c)
        try:
            A1_p_c = [norm(vmin=vmin, vmax=vmax)(p_c) for p_c in A1_p_c]
        except TypeError:
            # e.g. for MidNorm or mpl.colors.LogNorm, or if vmin, vmax already set
            A1_p_c = [norm(p_c) for p_c in A1_p_c]
        # Colour map etc
        if cmap == "fade":
            A1_colour[:, 3] = A1_p_c
        elif cmap == "tweak_luminosity":
            A1_colour = [tweak_luminosity(c, p_c) for p_c in A1_p_c]
        elif cmap is not None:
            A1_colour = cmap(A1_p_c)
    # Alpha
    A1_colour[:, 3] *= alpha

    # Linewidths
    if A1_lw is None:
        A1_lw = np.full(num_seg, lw)

    # Linestyle, done crudely by setting which segments to actually plot
    if ls == "-":
        A1_sel_seg = np.arange(num_seg)
    else:
        # Auto steps
        if ls == "--":
            ls_step_on = 0.03
            ls_step_off = 0.04
        elif ls == ":":
            ls_step_on = 0.01
            ls_step_off = 0.03

        # Projected distance between steps
        A1_d = np.cumsum(root_sum_sq([A1_x[:-1] - A1_x[1:], A1_y[:-1] - A1_y[1:]]))

        # Scale by the axis limits (unset limits are just 0, 1, so no effect)
        scale = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
        d_on = ls_step_on * scale
        d_off = ls_step_off * scale

        # Select subsets of ~equal-spaced segments as a crude workaround
        A1_sel_seg = np.where(A1_d % (d_on + d_off) < d_on)[0]

    # Plot
    lc = mpl.collections.LineCollection(
        segments=segments[A1_sel_seg],
        colors=A1_colour[A1_sel_seg],
        linewidths=A1_lw[A1_sel_seg],
        zorder=zorder,
        label=label,
    )
    ax.add_collection(lc)


def plot_orbit(
    o,
    A1_pos_0=np.zeros(3),
    unit=1,
    proj="xy",
    num_step=None,
    gradient=None,
    c="k",
    lw=1.7,
    alpha=1,
    zorder=None,
    label=None,
    ls="-",
    ls_step_on=None,
    ls_step_off=None,
    A1_ax_lim=None,
    mark_periapsis=False,
    plot_apse_line=False,
    eval_sel_rem=None,
    ax=None,
):
    """Plot an orbit.

    Parameters
    ----------
    o : Orbit
        An object containing orbital parameters. Required attributes: a, e, i,
        Omega, omega. Angles default to 0. Optional attributes: nu.

    A1_pos_0 : [float] (opt.)
        The x, y, z coordinates of the origin (the primary) (m).

    unit : float (opt.)
        The distance unit to plot in (m).

    proj : str (opt.)
        The plotting projection plane: xy, xz, zy, or yz.

    num_step : int (opt.)
        The number of steps to plot.

    gradient : str (opt.)
        "nu"        A gradient (of opacity) with true anomaly.
        "z"         A gradient (of linewidth) with z coordinate.

    c, lw, alpha, zorder, label : str, str, float, float, float, str (opt.)
        The colour, linewidth, opacity, zorder, and label.

    ls : str (opt.)
        The linestyle, done as a simple workaround, choose from: "-", "--", ":".
        Requires axis limits to be set or provided to be reliably effective.

    ls_step_on, ls_step_off : float, float (opt.)
        Instead of ls, set the linestyle directly with the lengths of the on and
        off segments, as fractions of the axis side.

    A1_ax_lim : [float] (opt.)
        The axis limits [x_min, x_max, y_min, y_max]. Allows ls options to work
        well, and avoids plotting outside the limits. Set to `"auto"` to use the
        figure values.

    mark_periapsis : bool (opt.)
        If True then plot a marker at the periapsis.

    plot_apse_line : bool (opt.)
        If True then plot the line connecting apoapsis and periapsis.

    eval_sel_rem : str (opt.)
        If provided, then evaluate this string as the argument for np.where(),
        to select part(s) of the line to remove, e.g. "A2_pos[2] < 0".

    ax : plt axes (opt.)
        The axes on which to plot. Defaults to the current axes.
    """
    # Defaults
    if ax is None:
        ax = plt.gca()
    if o.i is None:
        o.i = 0
    if o.Omega is None:
        o.Omega = 0
    if o.omega is None:
        o.omega = 0

    # Defaults
    if num_step is None:
        # Need finer steps for non-solid linestyles
        if ls == "-":
            num_step = 300
        elif o.e < 1:
            num_step = int(1e5)
        elif o.e >= 1:
            num_step = int(1e6)

    # List of true anomalies
    if o.is_parabolic:
        # Sample parabolic orbit only by true anomaly
        A1_nu = np.linspace(0, 2 * pi, num_step)
    else:
        # Want even number of steps to split between two sampling angles
        if num_step % 2 == 1:
            num_step += 1
        num_per_angle = int(num_step / 2)

        # Sample elliptical orbit by true and also eccentric anomaly
        if o.e < 1:
            A1_sample = np.linspace(0, 2 * pi, num_per_angle)

            # True anomaly and eccentric anomaly
            A1_nu = np.sort(np.append(A1_sample, [o.nu_from_E(E) for E in A1_sample]))
        # Sample hyperbolic orbit by true and also mean anomaly
        else:
            # Only angles between the aymptotes
            nu_asymp = np.arccos(-1.0 / o.e) * 0.995
            A1_sample = np.linspace(-nu_asymp, nu_asymp, num_per_angle)

            # True anomaly and mean anomaly
            A1_nu = np.sort(np.append(A1_sample, [o.nu_from_M(M) for M in A1_sample]))

    # Sort to start at the current position around the orbit, if known
    if o.nu is not None and not np.isnan(o.nu):
        A1_nu = np.append(A1_nu[A1_nu > o.nu], A1_nu[A1_nu <= o.nu])
        # Append current position
        if o.e < 1:
            A1_nu = np.append(o.nu, A1_nu)
            A1_nu = np.append(A1_nu, o.nu)
            num_step += 2

    # Calculate position at each true anomaly
    A2_pos = np.empty((num_step, 3))
    for i_step, nu in enumerate(A1_nu):
        A2_pos[i_step] = o.pos_vel_from_nu(nu, do_vel=False)

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
    else:
        raise Exception("Invalid proj", proj)

    # Relative to the primary # barycentre
    A2_pos[:, x_co] += A1_pos_0[x_co]  # + o.A1_pos_bary[x_co]
    A2_pos[:, y_co] += A1_pos_0[y_co]  # + o.A1_pos_bary[y_co]

    # Gradients
    num_seg = len(A2_pos)
    A1_p_c = None
    cmap = None
    A1_lw = None
    if gradient == "nu":
        # Projected distance between steps around orbit
        A1_d = (
            np.cumsum(
                root_sum_sq(
                    [
                        A2_pos[:-1, x_co] - A2_pos[1:, x_co],
                        A2_pos[:-1, y_co] - A2_pos[1:, y_co],
                    ]
                )
            )
            / unit
        )

        # Distance scale
        if o.e < 1:
            d_scale = np.amax(A1_d)
            # Not too large-scale for nearly-unbound orbits
            if o.e > 0.995:
                d_scale = 2 / (1 / d_scale + 1 / (400 * o.q / unit))
        else:
            d_scale = 200 * o.e * o.q / unit

        # Use distance around orbit instead of true anomaly for more even fading
        A1_p_c = np.minimum(0.02 + 2.1 * A1_d / d_scale, 1)
        cmap = "fade"
    elif gradient == "z":
        # cmap == "tweak_luminosity"
        # A1_z_scale = A2_pos[:, z_co] / abs(o.a) * 0.3 + 0.5
        A1_z_scale = A2_pos[:, z_co] / abs(o.a) * 3
        if np.amin(A1_z_scale) < -0.9:
            A1_z_scale += -0.9 - np.amin(A1_z_scale)
        A1_lw = np.full(num_seg, lw) * (1 + A1_z_scale)
    else:
        A1_lw = np.full(num_seg, lw)

    # Remove custom part(s) of the line
    if eval_sel_rem is not None:
        A1_sel_rem = np.where(eval(eval_sel_rem))[0]
        A2_pos[A1_sel_rem, :] = np.nan
        if A1_p_c is not None:
            A1_p_c[A1_sel_rem] = np.nan
        if A1_lw is not None:
            A1_lw[A1_sel_rem] = np.nan

    # Plot
    plot_lines(
        A1_x=A2_pos[:, x_co] / unit,
        A1_y=A2_pos[:, y_co] / unit,
        c=c,
        A1_p_c=A1_p_c,
        cmap=cmap,
        A1_lw=A1_lw,
        ls=ls,
        ls_step_on=ls_step_on,
        ls_step_off=ls_step_off,
        A1_ax_lim=A1_ax_lim,
        alpha=alpha,
        zorder=zorder,
        label=label,
        ax=ax,
    )

    # Annotations
    if mark_periapsis:
        A1_pos_q = A1_pos_0 + o.A1_pos_q
        for marker in ["1", "2"]:
            ax.scatter(
                [A1_pos_q[x_co] / unit],
                [A1_pos_q[y_co] / unit],
                marker=marker,
                s=4**2,
                c=[c],
                edgecolor="none",
                lw=lw * 0.7,
                alpha=alpha * 0.7,
            )
    if plot_apse_line:
        A1_pos_q = A1_pos_0 + o.A1_pos_q
        A1_pos_Q = A1_pos_0 + o.A1_pos_Q
        ax.plot(
            [A1_pos_q[x_co] / unit, A1_pos_Q[x_co] / unit],
            [A1_pos_q[y_co] / unit, A1_pos_Q[y_co] / unit],
            c=c,
            lw=lw * 0.7,
            alpha=alpha * 0.7,
        )


def plot_hist(
    ax,
    A1_hist,
    A1_bin_edge,
    c="k",
    ls="-",
    lw=1.7,
    alpha=1,
    zorder=None,
    label=None,
    smooth=0,
):
    """Plot a nice histogram with no vertical lines between bins.

    Parameters
    ----------
    ax : Axes
        The plot axes.

    A1_hist, A1_bin_edge : [float]
        The histogram and bin arrays returned by numpy.histogram(): a list of
        number counts or density values in each bin and a list of bin edges
        including both outside values.

    c, ls, lw, alpha, zorder : str, str, str, float, float (opt.)
        The colour, linestyle, linewidth, opacity, and zorder. Default solid
        black line.

    label : str (opt.)
        The label, if required.

    smooth : float (opt.)
        If >0, then do a crude "smoothing" between the bins for sometimes-
        nicer plotting. The value (<0.5) sets how much of the bin to smooth.
        Setting smooth=1 plots a normal line through bin centres.
    """
    if smooth is None:
        smooth = 0

    # Plot a normal line through the bin centres
    if smooth == 1:
        A1_bin_mid = (A1_bin_edge[1:] + A1_bin_edge[:-1]) / 2
        ax.plot(
            A1_bin_mid,
            A1_hist,
            c=c,
            ls=ls,
            lw=lw,
            alpha=alpha,
            zorder=zorder,
            label=label,
        )
        return

    # Append values to complete the plot with vertical lines at both ends
    A1_hist = np.append(-np.inf, A1_hist)
    A1_hist = np.append(A1_hist, -np.inf)
    A1_bin_edge = np.append(A1_bin_edge, A1_bin_edge[-1])

    # Plot
    if smooth > 0:
        # Crudely smooth the step from bin to bin
        # Bin widths
        A1_bin_width = A1_bin_edge[1:] - A1_bin_edge[:-1]
        A1_bin_width = np.append(A1_bin_width, A1_bin_width[-1])

        # Extra values for normal stepped histogram first
        A1_bin_edge = np.repeat(A1_bin_edge, 2)[:-1]
        A1_hist = np.repeat(A1_hist, 2)[1:]

        # Tweaked values for smoothing
        A1_bin_edge[::2] = A1_bin_edge[::2] - smooth * A1_bin_width
        A1_bin_edge[1::2] = A1_bin_edge[1::2] + smooth * A1_bin_width[1:]

        # # Smoothed spline
        # A1_x_smooth = np.linspace(A1_bin_edge[0], A1_bin_edge[1], 500)
        # spl = make_interp_spline(A1_bin_edge, A1_hist, k=3)
        # A1_y_smooth = spl(A1_x_smooth)

        # Plot
        ax.plot(
            A1_bin_edge,
            A1_hist,
            c=c,
            ls=ls,
            lw=lw,
            alpha=alpha,
            zorder=zorder,
            label=label,
        )
    else:
        # Normal stepped histogram
        ax.plot(
            A1_bin_edge,
            A1_hist,
            c=c,
            ls=ls,
            lw=lw,
            drawstyle="steps",
            alpha=alpha,
            zorder=zorder,
            label=label,
        )


def plot_gradient_line(
    A1_x,
    A1_y,
    A1_c,
    cmap,
    vmin=None,
    vmax=None,
    norm=mpl.colors.Normalize,
    x_scale=None,
    y_scale=None,
    lw=1.7,
    alpha=1,
    zorder=None,
    ax=None,
):
    """Plot a line with a colour gradient following values with a colour map.

    Parameters
    ----------
    A1_x, A1_y : [float]
        The x and y coordinates to plot.

    A1_c : [float]
        The colour values.

    cmap  : Colormap
        The colour map to use.

    vmin, vmax, norm : float, float, Normalize (opt.)
        The min and max values and the normalisation for the colour map.

    x_scale, y_scale : str (opt.)
        Set to "log" if the plotting axis will be logarithmic.

    lw, alpha, zorder : float, float, float (opt.)
        The linewidth, opacity, and zorder.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.
    """
    if ax is None:
        ax = plt.gca()

    # Normalise the colour values to between 0 and 1 for the colour map
    if vmin is None:
        vmin = np.amin(A1_c)
    if vmax is None:
        vmax = np.amax(A1_c)
    try:
        A1_c_norm = [norm(vmin=vmin, vmax=vmax)(c) for c in A1_c]
    except TypeError:
        # e.g. for MidNorm or mpl.colors.LogNorm, or if vmin, vmax already set
        A1_c_norm = [norm(c) for c in A1_c]

    # Interpolate, depending on how many points are already in the line
    n = len(A1_x)
    if n < 10:
        n_intp = 200
    elif n < 100:
        n_intp = 20
    else:
        n_intp = 2

    # Plot between one pair of values at a time
    for i in range(n - 1):
        # Interpolate
        if x_scale == "log":
            A1_x_intp = np.geomspace(A1_x[i], A1_x[i + 1], n_intp + 1)
        else:
            A1_x_intp = np.linspace(A1_x[i], A1_x[i + 1], n_intp + 1)
        if y_scale == "log":
            A1_y_intp = np.geomspace(A1_y[i], A1_y[i + 1], n_intp + 1)
        else:
            A1_y_intp = np.linspace(A1_y[i], A1_y[i + 1], n_intp + 1)
        A1_c_intp = np.linspace(A1_c_norm[i], A1_c_norm[i + 1], n_intp + 1)

        # Set colour for each mini line
        ax.set_prop_cycle(color=[cmap(c) for c in A1_c_intp])

        # Plot each mini segment
        for j in range(n_intp):
            ax.plot(
                [A1_x_intp[j], A1_x_intp[j + 1]],
                [A1_y_intp[j], A1_y_intp[j + 1]],
                lw=lw,
                alpha=alpha,
                zorder=zorder,
            )


def plot_opacity_gradient_line(
    A1_x,
    A1_y,
    alpha_0,
    alpha_1,
    n_intp=None,
    x_scale=None,
    y_scale=None,
    c="k",
    lw=1.7,
    marker=None,
    zorder=None,
    ax=None,
):
    """Plot a line with an opacity gradient.

    Parameters
    ----------
    A1_x, A1_y : [float]
        The x and y coordinates to plot.

    alpha_0, alpha_1 : float (opt.)
        The starting and final opacity.

    n_intp : int (opt.)
        The number of interpolation steps to do between each data point.

    x_scale, y_scale : str (opt.)
        Set to "log" if the plotting axis will be logarithmic.

    c, lw, zorder : str, float, float (opt.)
        The colour, linewidth, and zorder.

    marker : str (opt.)
        If provided, then also plot markers with the appropriate opacities.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.
    """
    if ax is None:
        ax = plt.gca()

    # Interpolate, depending on how many points are already in the line
    n = len(A1_x)
    if n_intp is None:
        if n < 10:
            n_intp = 200
        elif n < 100:
            n_intp = 20
        else:
            n_intp = 2

    # Opacities
    A1_alpha_intp = np.linspace(alpha_0, alpha_1, n * n_intp)

    # Plot between one pair of values at a time
    for i in range(n - 1):
        # Interpolate
        if x_scale == "log":
            A1_x_intp = np.geomspace(A1_x[i], A1_x[i + 1], n_intp + 1)
        else:
            A1_x_intp = np.linspace(A1_x[i], A1_x[i + 1], n_intp + 1)
        if y_scale == "log":
            A1_y_intp = np.geomspace(A1_y[i], A1_y[i + 1], n_intp + 1)
        else:
            A1_y_intp = np.linspace(A1_y[i], A1_y[i + 1], n_intp + 1)

        # Plot each mini segment
        for j in range(n_intp):
            ax.plot(
                [A1_x_intp[j], A1_x_intp[j + 1]],
                [A1_y_intp[j], A1_y_intp[j + 1]],
                c=c,
                lw=lw,
                alpha=A1_alpha_intp[i * n_intp + j],
                zorder=zorder,
            )

    # Plot the markers
    if marker is not None:
        A1_alpha = np.linspace(alpha_0, alpha_1, n)

        for i in range(n):
            ax.scatter(
                A1_x[i], A1_y[i], marker=marker, c=c, alpha=A1_alpha[i], zorder=zorder
            )


def plot_angle(
    A1_x, A1_y, offset=1, c="k", ls="-", lw=1.7, alpha=1, zorder=None, ax=None
):
    """Plot an angle arc marking between two lines.

    Parameters
    ----------
    A1_x, A1_y : [float, float, float]
        The coordinates of the three points, with the pivot point in the middle.

    offset : float (opt.)
        The size of the arc.

    c, ls, lw, alpha, zorder : str, str, float, float, float (opt.)
        The colour, linestyle, linewidth, opacity, and zorder.
        Defaults to a solid black line.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.
    """
    if ax is None:
        ax = plt.gca()

    pivot = (A1_x[1], A1_y[1])

    # Angles wrt x axis
    angle_1 = np.arctan2(A1_y[0] - pivot[1], A1_x[0] - pivot[0])
    angle_2 = np.arctan2(A1_y[2] - pivot[1], A1_x[2] - pivot[0])

    arc = Arc(
        pivot,
        offset,
        offset,
        angle=0,
        theta1=min(angle_1, angle_2) * rad_to_deg,
        theta2=max(angle_1, angle_2) * rad_to_deg,
        color=c,
        ls=ls,
        lw=lw,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(arc)


def single_arrow(
    xy_1, xy_2, ax=None, c="k", ls="-", lw=1.7, alpha=1, head_width=5, head_length=8
):
    """Plot a single-headed arrow between two points.

    Parameters
    ----------
    xy_1, xy_2 : (float, float)
        The x, y coordinates of the two points.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.

    c, ls, lw, alpha : str (opt.)
        The colour, linestyle, linewidth, and opacity. Default solid black line.
    """
    if ax is None:
        ax = plt.gca()

    arrow = FancyArrowPatch(
        xy_1,
        xy_2,
        arrowstyle="->, head_width=%f, head_length=%f" % (head_width, head_length),
        shrinkA=0,
        shrinkB=0,
        color=c,
        ls=ls,
        lw=lw,
        alpha=alpha,
    )
    ax.add_patch(arrow)


def double_arrow(xy_1, xy_2, ax=None, c="k", ls="-", lw=1.7, alpha=1):
    """Plot a double-headed arrow between two points.

    Parameters
    ----------
    xy_1, xy_2 : (float, float)
        The x, y coordinates of the two points.

    ax : Axes (opt.)
        The axes on which to plot. Defaults to the current axes.

    c, ls, lw, alpha : str (opt.)
        The colour, linestyle, linewidth, and opacity. Default solid black line.
    """
    if ax is None:
        ax = plt.gca()

    arrow = FancyArrowPatch(
        xy_1,
        xy_2,
        arrowstyle="<->, head_width=5, head_length=8",
        shrinkA=0,
        shrinkB=0,
        color=c,
        ls=ls,
        lw=lw,
        alpha=alpha,
    )
    ax.add_patch(arrow)


def text_outlined(x, y, text, outl_w=6, outl_c="w", ax=None, **kwargs):
    """Standard text, with an outline border."""
    if ax is None:
        ax = plt.gca()

    txt = ax.text(x, y, text, **kwargs)

    if outl_w > 0:
        txt.set_path_effects(
            [path_effects.withStroke(linewidth=outl_w, foreground=outl_c)]
        )


def add_text_loc(
    text, loc, c="k", pad=0.03, outl_w=6, outl_c="w", fontsize=24, xy_ratio=1, ax=None
):
    """Add text to a plot.

    Parameters
    ----------
    text : str
        The text to add.

    loc : str
        The location for the text: "upper/lower/center left/right/center".

    c : str (opt.)
        The text colour.

    pad : float (opt.)
        The distance of the text from the axis edges, as a fraction of the side.

    outl_w, outl_c : float, str (opt.)
        The width and colour of a background outline.

    fontsize : float (opt.)
        The font size.

    xy_ratio : float (opt.)
        The ratio of the size of the x and y axis, used to scale the text
        position for non-square axes (= x_side / y_side).

    ax : Axes (opt.)
        The plot axes.
    """
    if ax is None:
        ax = plt.gca()

    # Location
    if "upper" in loc:
        y = 1 - pad
        va = "top"
    elif "lower" in loc:
        y = pad
        va = "bottom"
    elif loc[:6] == "center":
        y = 0.5
        va = "center"
    else:
        raise Exception("Invalid loc ", loc)

    if "left" in loc:
        x = pad
        ha = "left"
    elif "right" in loc:
        x = 1 - pad
        ha = "right"
    elif loc[-6:] == "center":
        x = 0.5
        ha = "center"
    else:
        raise Exception("Invalid loc ", loc)

    # Adjust for unequal axes
    if xy_ratio < 1:
        if x < 0.5:
            x /= xy_ratio
        elif x > 0.5:
            x = 1 - (1 - x) / xy_ratio
    elif xy_ratio > 1:
        if y < 0.5:
            y *= xy_ratio
        elif y > 0.5:
            y = 1 - (1 - y) * xy_ratio

    # Add text
    text_outlined(
        x,
        y,
        text,
        outl_w=outl_w,
        outl_c=outl_c,
        ax=ax,
        transform=ax.transAxes,
        fontsize=fontsize,
        color=c,
        ha=ha,
        va=va,
    )


# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  jkeger.py  ====\n")
