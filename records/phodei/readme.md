Project `phodei`
================

## Origin of Mars’s moons by disruptive partial capture of an asteroid


Files
-----
+ `Ma_xp_A2000_n65_r11_v00.yml` Base input parameter file for SWIFT simulations.
+ `output_list_120ks.txt` The list of times for snapshots to be created.
+ `select_output.yml` The particle data fields to be output in snapshots.
+ `rebound/`
    + `Makefile` Makefile for compiling the REBOUND script.
    + `simulation.c` Script for running the REBOUND simulations.


SWIFT
-----
+ SHA b153f18f800051e07e5e3ae204e644dc43b5161b, subtask_speedup
+ `./configure --with-hydro=planetary --with-equation-of-state=planetary --with-gravity=basic --with-ext-potential=point-mass --with-tbbmalloc --with-parmetis --enable-ipo`
+ Varying input parameters for different simulations:
    + `file_name`, `basename`, etc.
    + `h_max` and `max_physical_baryon_softening` change with particle mass (i.e. with resolution and/or parent asteroid mass) as detailed in `objects_phodei.py`.
