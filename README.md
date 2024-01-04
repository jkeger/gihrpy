Giant Impacts at High Resolution (GIHR) Python
==============================================

Jacob Kegerreis  jacob.kegerreis@durham.ac.uk

Python code for GIHR simulation projects and related bits and pieces. Mostly
things like making initial conditions or analysing and plotting results.
See also [SWIFT](swiftsim.com) and [WoMa](https://github.com/srbonilla/WoMa).

This code is provided in the hopes that it may be useful and to improve the
openness and reproducibility of published projects. However, please note that
this code has grown organically and often messily over many years and many
different projects, and was not designed primarily with other users in mind. As
such, while there is hopefully sufficient documentation for the code to be
usable and intelligible, there is a lot that could be refactored and tidied.


Files
-----
+ `gihr.py` Main script, see its doc string for usage.

+ `gihrpy/`
    + `init_cond.py` Functions for generating and plotting initial conditions.
    + `snapshot_plots.py` Functions for plotting standard snapshot data, etc.
    + `scenario_plots.py` Functions for plotting accumulated scenario results, etc.
    + `rebound_plots.py` Functions for plotting rebound results.
    + `custom_plots.py` Functions for customisable plots like tiled figures, etc.

    + `objects.py` Accumulated input-info objects from all projects for easy access.
    + `objects_*.py` Input objects for specific projects.

    + `classes.py` Classes for plotting and initial conditions input objects.
    + `classes_sim.py` Classes for main simulation input objects.
    + `classes_reb.py` Classes for rebound input objects.

    + `utilities.py` General utility functions for e.g. loading data and plotting extras.
    + `utilities_reb.py` Rebound utility functions for e.g. loading data.
    + `utilities_*.py` Extra utilities for specific projects, like custom plotting tweaks.

+ `orbits.py` Stand-alone module for computing and plotting Keplerian orbits.
+ `jkeger.py` Various generic utilities and constants, etc.
+ `tables/`
    + `cmap_*.pkl` Generated colour maps.

+ `records/`
    + `*/` Extra recorded info for specific projects, e.g. simulation input parameters.
        + `readme.md` Details about this specific project and any other record files.

+ `README.md` This file.
+ `LICENSE.txt` GNU general public license v3+.
+ `format.sh` Simple formatting script.


Notation etc.
-------------
+ PEP8 is followed in most cases, auto-formatted with [black](https://github.com/psf/black).
+ Arrays are explicitly labelled with a prefix `A1_`, or `An_` for an `n`-dimensional array.
+ Dictionaries are labelled with a prefix `Di_`.
+ File paths are labelled with a prefix `Fp_`.
+ Units are SI unless specified.
