/**
GIHR REBOUND simulation source file.

Contains various utilities and modified functions for e.g. units and I/O.

Parameters
----------
sim : str
    The name of the simulation to run. There should exist a folder of this name
    with the initial conditions and config files (see load functions) inside it.

-h, --help
    Print help information and exit.

-c <str>, --config=<str>
    The name of the config file to load, assumed to be inside the simulation
    directory (or a relative path from there). Default `config.txt`.

-r, --restart
    Restart the simulation from a previously made restart file, `restart.bin`.

-n <int>, --n_per_split=<int>
-i <int>, --i_split=<int>
    For splitting up test particles between runs for crude parallelisation, the
    number of test particles per run and the index of this run. e.g. `-n 5 -i 0`
    with 25 test particles: simulate only the first 5 test particles, then
    `-n 5 -i 1`: simulate the next 5, up to `-n -i 4`: simulate the last 5.

-C, --clean
    Remove any existing output files for a clean start.

-d, --dry
    Do a dry run without writing any output files.

-s <int>, --steps=<int>
    Run for only this many steps.

-w <float>, --wall_max_h=<float>
    The wall-clock maximum runtime in hours, defaults to 71.5 h.

-t, --test
    Run temporary testing code.
**/
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rebound.h"

// Units (for internal G = 1)
const double m_unit = 1.988475e30;   // kg (M_Sun)
const double l_unit = 1.4959787e11;  // m (AU)
const double t_unit = 5.02254803e6;  // s (yr/2pi)
const double v_unit = 2.97852542e4;  // m/s
const double year = 3.15576e7;       // s

// Simulation name and input files
char sim[999];
char file_config[999];
char file_init_cond[999];
char file_out[999];
char file_restart[999];
char file_restart_prev[999];
char file_collisions[999];

// Split test particles between runs
int num_test;
int n_per_split = 0;
int i_split = 0;

// Output times
double t_start;
double t_end;
double t_out_step;
double t_end_rel;
double t_out_step_1;
double t_2;
double t_out_step_2;
double t_3;
double t_out_step_3;

// Restarts
int do_restart = 0;
int restart_step = 1;
double t_wall_restart;

// Relative boundary removal
int id_rel_prim;
double rel_distance_max;

// Oblate potential
int oblate_id;
double oblate_J2;
double oblate_obliquity;

// Collisions
int collision_mode;
double coeff_restitution;
double t_collision_delay;

// Reference frame-shift
int id_frame_shift;

// Misc
int do_clean = 0;
int do_dry_run = 0;
int steps_end = -1;
float wall_max_h = 71.5;
int is_first_heartbeat = 1;
int is_first_restart = 1;
int is_first_collision = 1;
int i_output_start = 0;

/*
    Print an error message and exit.
*/
#define error(message, ...)                           \
    {                                                 \
        fflush(stdout);                               \
        fprintf(stderr, "Error (%s): \n  ", sim);     \
        fprintf(stderr, message "\n", ##__VA_ARGS__); \
        exit(1);                                      \
    }

/*
    Print a warning message.
*/
#define warning(message, ...)                         \
    {                                                 \
        fflush(stdout);                               \
        fprintf(stderr, "Warning (%s): \n  ", sim);   \
        fprintf(stderr, message "\n", ##__VA_ARGS__); \
    }

/*
    Skip a line while reading a file.
*/
int skip_line(FILE* f) {
    int c;

    // Read each character until reaching the end of the line or file
    do {
        c = fgetc(f);
    } while ((c != '\n') && (c != EOF));

    return c;
}

/*
    Skip n lines while reading a file.
*/
int skip_lines(FILE* f, int n) {
    int c;

    for (int i = 0; i < n; i++) c = skip_line(f);

    return c;
}

/*
    Load configuration parameters from a text file.

    Set `-1` for default values or to skip an optional parameter.

    File contents
    -------------
    # HEADER (8 lines)
    num_picle  num_test                                     (-  -)
    file_init_cond  file_out_stem  file_restart             (-  -  -)
    t_start  t_end  t_out_step  t_wall_restart              (yr  yr  yr  h)
    t_2  t_out_step_2  t_3  t_out_step_3                    (yr  yr  yr  yr)
    id_rel_prim  rel_distance_max                           (-  m)
    oblate_id  oblate_J2  oblate_obliquity                  (-  -  rad)
    collision_mode  coeff_restitution  t_collision_delay    (-  -  yr)
    id_frame_shift                                          (-)

    Contents are as defined below. The provided file paths are assumed to be
    relative to the simulation directory (like --config).

    Parameters
    ----------
    filename : str
        The path to the file to load.

    Returns
    -------
    num_picle : int
        The number of particles.

    Sets
    ----
    num_test : int
        The number of test particles (or mark test particles with zero mass).
        Assumes non-test particles listed first, override the rest to zero mass.

    file_init_cond : str
        The initial conditions file path.

    file_out_stem : str
        The output file path, not including the .txt extension.

    file_restart : str
        The restart file path.

    t_start : double
        The start time (s).

    t_end : double
        The end time (s).

    t_out_step : double
        The time between outputs (s).

    t_wall_restart : double
        The wall-clock time between writing restart files (h).

    t_end_rel : double
        The integration end time, relative to the start time (s).

    t_2(_3) : double
        The time to change to the second (third) output step (s).

    t_out_step_2(_3) : double
        The second (third) time between outputs, from t_2 (t_3) onwards (s).

    id_rel_prim : int
        The ID of the designated primary particle to remove orbiting particles
        that get too close or too far from it.

    rel_distance_max : double
        The maximum distance from the designated primary particle beyond which
        particles will be removed.

    oblate_id : int
        The ID of the oblate particle.

    oblate_J2 : double
        The J2 moment of the oblate particle.

    oblate_obliquity : double
        The obliquity (rad) of the oblate particle.

    collision_mode : int
        The type of inter-particle collisions to use: -1 = no collisions;
        0 = print collision info only; 1 = hard-sphere collisions; 2 = merge.

    coeff_restitution : double
        The coefficient of restitution for hard-sphere collisions.

    t_collision_delay : double
        Delay collisions to only be allowed after this time has elapsed (s),
        e.g. to allow particles to separate from close initial conditions.

    id_frame_shift : int
        The ID of a particle to which the reference frame should be shifted
        every heartbeat. Can yield a big speed-up for e.g. a planet with
        orbiting particles that are together all orbiting a star.
*/
int load_config_parameters(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) error("Failed to open file '%s'", filename);

    // Header
    int c = skip_lines(f, 8);

    // Data
    int num_picle;
    c = fscanf(f, "%d %d", &num_picle, &num_test);
    if (c != 2) error("Failed to read num_picle, num_test '%s' (returned %d)", filename, c);

    char file_init_cond_in[999];
    char file_out_stem_in[999];
    char file_restart_in[999];
    c = fscanf(f, "%s %s %s", file_init_cond_in, file_out_stem_in, file_restart_in);
    if (c != 3)
        error(
            "Failed to read file_init_cond, file_out_stem, file_restart '%s' (returned %d)",
            filename, c);
    // Set relative to the simulation directory
    sprintf(file_init_cond, "%s/%s", sim, file_init_cond_in);
    if (n_per_split != 0) {
        // Label by split number
        sprintf(file_out, "%s/%s_%d.txt", sim, file_out_stem_in, i_split);
    } else {
        sprintf(file_out, "%s/%s.txt", sim, file_out_stem_in);
    }
    sprintf(file_restart, "%s/%s", sim, file_restart_in);
    sprintf(file_restart_prev, "%s.prev", file_restart);

    c = fscanf(f, "%lf %lf %lf %lf", &t_start, &t_end, &t_out_step, &t_wall_restart);
    if (c != 4)
        error(
            "Failed to read t_start, t_end, t_out_step, t_wall_restart '%s' (returned %d)",
            filename, c);

    c = fscanf(f, "%lf %lf %lf %lf", &t_2, &t_out_step_2, &t_3, &t_out_step_3);
    if (c != 4)
        error(
            "Failed to read t_2, t_out_step_2, t_3, t_out_step_3 '%s' (returned %d)", filename, c);

    c = fscanf(f, "%d %lf", &id_rel_prim, &rel_distance_max);
    if (c != 2)
        error("Failed to read id_rel_prim, rel_distance_max '%s' (returned %d)", filename, c);

    c = fscanf(f, "%d %lf %lf", &oblate_id, &oblate_J2, &oblate_obliquity);
    if (c != 3) error("Failed to read J2 parameters '%s' (returned %d)", filename, c);

    c = fscanf(f, "%d %lf %lf", &collision_mode, &coeff_restitution, &t_collision_delay);
    if (c != 3) error("Failed to read collision parameters '%s' (returned %d)", filename, c);

    c = fscanf(f, "%d", &id_frame_shift);
    if (c != 1) error("Failed to read id_frame_shift '%s' (returned %d)", filename, c);

    fclose(f);

    // Convert to seconds
    t_start *= year;
    t_end *= year;
    t_out_step *= year;
    t_wall_restart *= 3600.0;
    t_out_step_1 = t_out_step;
    t_2 *= year;
    t_out_step_2 *= year;
    t_3 *= year;
    t_out_step_3 *= year;
    t_collision_delay *= year;

    // Relative integration time
    t_end_rel = t_end - t_start;

    // Collisions output file
    sprintf(file_collisions, "%s/collisions.txt", sim);

    // Clean old files
    if (do_clean) {
        // Back up the cleaned files just in case
        char file_old[999];
        sprintf(file_old, "%s.old", file_out);
        rename(file_out, file_old);
        sprintf(file_old, "%s.old", file_restart);
        rename(file_restart, file_old);
        sprintf(file_old, "%s.old", file_restart_prev);
        rename(file_restart_prev, file_old);
        sprintf(file_old, "%s.old", file_collisions);
        rename(file_collisions, file_old);

        // Also any previous output files
        char file_remove[999];
        sprintf(file_remove, "%s.prev", file_out);
        sprintf(file_old, "%s.old", file_remove);
        rename(file_remove, file_old);
        sprintf(file_remove, "%s.prev", file_collisions);
        sprintf(file_old, "%s.old", file_remove);
        rename(file_remove, file_old);
    }

    return num_picle;
}

/*
    Load the initial conditions from a text file.

    If num_test > 0, then assume all non-test particles are listed first. And/or
    test particles can be set by zero mass.

    File contents
    -------------
    # HEADER (2 lines)
    num_picle                                       (-)
    m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0    (kg  m  m  m  m  m/s  m/s  m/s)
    m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                 ...                  ...
    m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n

    Parameters
    ----------
    filename : str
        The path to the file to load.

    num_picle : int
        The number of particles.

    A1_m : [double]
        The particle masses (kg).

    A1_R : [double]
        The particle radii (m).

    A2_pos : [[double]]
        The particle positions (m).

    A2_vel : [[double]]
        The particle velocities (m/s).
*/
void load_initial_conditions(
    const char* filename, int num_picle, double A1_m[num_picle], double A1_R[num_picle],
    double A2_pos[num_picle][3], double A2_vel[num_picle][3]) {
    FILE* f = fopen(filename, "r");
    if (!f) error("Failed to open file '%s'", filename);

    // Header
    int c = skip_lines(f, 2);

    // Data
    c = fscanf(f, "%d", &num_picle);
    if (c != 1) error("Failed to read num_picle '%s' (returned %d)", filename, c);

    // Load particle data
    for (int i = 0; i < num_picle; i++) {
        c = fscanf(
            f, "%lf %lf %lf %lf %lf %lf %lf %lf", &A1_m[i], &A1_R[i], &A2_pos[i][0], &A2_pos[i][1],
            &A2_pos[i][2], &A2_vel[i][0], &A2_vel[i][1], &A2_vel[i][2]);
        if (c != 8) error("Failed to read particle %d's data '%s' (returned %d)", i, filename, c);
    }

    fclose(f);
}

/*
    Resolve collisions, just write the info to a text file.

    File contents
    -------------
    # HEADER (1 line)
    t  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1  <same for _2>  r_prim  v_rel  (yr  SI...)
*/
int collision_resolve_write_info(struct reb_simulation* const r, struct reb_collision c) {
    if (do_dry_run) return 0;

    struct reb_particle p1 = r->particles[c.p1];
    struct reb_particle p2 = r->particles[c.p2];

    double time = (t_start + r->t * t_unit) / year;
    // Relative to the origin or the designated primary
    double x_0 = 0.0;
    double y_0 = 0.0;
    double z_0 = 0.0;
    double vx_0 = 0.0;
    double vy_0 = 0.0;
    double vz_0 = 0.0;
    if (r->id_rel_prim >= 0) {
        x_0 = r->particles[r->id_rel_prim].x;
        y_0 = r->particles[r->id_rel_prim].y;
        z_0 = r->particles[r->id_rel_prim].z;
        vx_0 = r->particles[r->id_rel_prim].vx;
        vy_0 = r->particles[r->id_rel_prim].vy;
        vz_0 = r->particles[r->id_rel_prim].vz;
    }
    double r_prim = sqrt(
        ((p1.x + p2.x) / 2 - x_0) * ((p1.x + p2.x) / 2 - x_0) +
        ((p1.y + p2.y) / 2 - y_0) * ((p1.y + p2.y) / 2 - y_0) +
        ((p1.z + p2.z) / 2 - z_0) * ((p1.z + p2.z) / 2 - z_0));
    // Relative velocity
    double v_rel = sqrt(
        (p2.vx - p1.vx) * (p2.vx - p1.vx) + (p2.vy - p1.vy) * (p2.vy - p1.vy) +
        (p2.vz - p1.vz) * (p2.vz - p1.vz));

    FILE* f = fopen(file_collisions, "a");
    if (!f) error("Failed to open file '%s'", file_collisions);

    // Header
    if ((is_first_collision) & (!do_restart)) {
        fprintf(
            f,
            "# t  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1  "
            "m_2  R_2  x_2  y_2  z_2  v_x_2  v_y_2  v_z_2  "
            "r_prim  v_rel  "
            "(yr  kg  m  m  m  m/s  m/s  m/s  kg  m  m  m  m/s  m/s  m/s  m  m/s)  "
            "# Pos & vel wrt primary \n");
    }

    // Append collision data
    fprintf(f, "%e ", time);
    fprintf(
        f, "%e %e %e %e %e %e %e %e ", p1.m * m_unit, p1.r * l_unit, (p1.x - x_0) * l_unit,
        (p1.y - y_0) * l_unit, (p1.z - z_0) * l_unit, (p1.vx - vx_0) * v_unit,
        (p1.vy - vy_0) * v_unit, (p1.vz - vz_0) * v_unit);
    fprintf(
        f, "%e %e %e %e %e %e %e %e ", p2.m * m_unit, p2.r * l_unit, (p2.x - x_0) * l_unit,
        (p2.y - y_0) * l_unit, (p2.z - z_0) * l_unit, (p2.vx - vx_0) * v_unit,
        (p2.vy - vy_0) * v_unit, (p2.vz - vz_0) * v_unit);
    fprintf(f, "%e %e \n", r_prim * l_unit, v_rel * v_unit);

    fclose(f);

    // Mark first collision done
    if (is_first_collision) is_first_collision = 0;

    // Don't remove either particle
    return 0;
}

/*
    Return the coefficient of restitution for resolving hard-sphere collisions,
    using the global-set fixed value.
*/
double coeff_restitution_fixed(const struct reb_simulation* const r, double v) {
    return coeff_restitution;
}

/*
    Resolve collisions, interact as hard spheres and write info.
*/
int collision_resolve_hardsphere_and_write(
    struct reb_simulation* const r, struct reb_collision c) {
    collision_resolve_write_info(r, c);
    int ret = reb_collision_resolve_hardsphere(r, c);
    return ret;
}

/*
    Resolve collisions, merge as hard spheres and write info.
*/
int collision_resolve_merge_and_write(struct reb_simulation* const r, struct reb_collision c) {
    collision_resolve_write_info(r, c);
    int ret = reb_collision_resolve_merge(r, c);
    return ret;
}

/*
    Checkpoint printing.

    Tweaked reb_output_timing(). Add steps, convert time units to years, etc.
*/
void reb_output_timing_v2(struct reb_simulation* r, const double tmax) {
    // Header
    if (is_first_heartbeat) {
        printf(
            "Steps      Particles   Time (yr) (frac)   Timestep (yr) CPU time (s) (tot, min)\n");
    }

    const int N = r->N;
#ifdef MPI
    int N_tot = 0;
    MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (r->mpi_id != 0) return;
#else
    int N_tot = N;
#endif
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double temp = tim.tv_sec + (tim.tv_usec / 1000000.0);
    if (r->output_timing_last == -1) {
        r->output_timing_last = temp;
    } else {
        printf("\r");
#ifdef PROFILING
        fputs("\033[A\033[2K", stdout);
        for (int i = 0; i <= PROFILING_CAT_NUM; i++) {
            fputs("\033[A\033[2K", stdout);
        }
#endif  // PROFILING
    }

    if (r->steps_done < 1e8) {
        printf("%-9llu", r->steps_done);
    } else {
        printf("%.3e", (double)r->steps_done);
    }
    printf("  N=%-8d", N_tot);
    if (r->integrator == REB_INTEGRATOR_SEI) {
        printf("  t=%-7f [orb]", (t_start / t_unit + r->t) * r->ri_sei.OMEGA / 2.0 / M_PI);
    } else {
        printf("  t=%-6f", (t_start + r->t * t_unit) / year);
    }
    if (tmax > 0) {
        printf(" (%4.1f%%)", r->t / tmax * 100.0);
    }
    printf("  dt=%-7f", r->dt * t_unit / year);
    if (r->walltime / 60.0 > 120) {
        printf("  wall=%.2f (%.1fh)", temp - r->output_timing_last, r->walltime / 3600.0);
    } else {
        printf("  wall=%.3f (%.3g)", temp - r->output_timing_last, r->walltime / 60.0);
    }
#ifdef PROFILING
    if (profiling_timing_initial == 0) {
        struct timeval tim;
        gettimeofday(&tim, NULL);
        profiling_timing_initial = tim.tv_sec + (tim.tv_usec / 1000000.0);
    }
    printf("\nCATEGORY       TIME \n");
    double _sum = 0;
    for (int i = 0; i <= PROFILING_CAT_NUM; i++) {
        switch (i) {
            case PROFILING_CAT_INTEGRATOR:
                printf("Integrator     ");
                break;
            case PROFILING_CAT_BOUNDARY:
                printf("Boundary check ");
                break;
            case PROFILING_CAT_GRAVITY:
                printf("Gravity/Forces ");
                break;
            case PROFILING_CAT_COLLISION:
                printf("Collisions     ");
                break;
#ifdef OPENGL
            case PROFILING_CAT_VISUALIZATION:
                printf("Visualization  ");
                break;
#endif  // OPENGL
            case PROFILING_CAT_NUM:
                printf("Other          ");
                break;
        }
        if (i == PROFILING_CAT_NUM) {
            printf(
                "%5.2f%%", (1. - _sum / (profiling_time_final - profiling_timing_initial)) * 100.);
        } else {
            printf(
                "%5.2f%%\n",
                profiling_time_sum[i] / (profiling_time_final - profiling_timing_initial) * 100.);
            _sum += profiling_time_sum[i];
        }
    }
#endif  // PROFILING
    printf("\n");
    fflush(stdout);
    r->output_timing_last = temp;
}

/*
    Write/append the output data to a text file.

    File contents
    -------------
    # t  m  R  x  y  z  v_x  v_y  v_z    (yr  kg  m  m  m  m  m/s  m/s  m/s)
    t_0  m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0
    t_0  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                    ...                    ...
    t_0  m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n
    t_1  m_0  R_0  x_0  y_0  z_0  v_x_0  v_y_0  v_z_0
    t_1  m_1  R_1  x_1  y_1  z_1  v_x_1  v_y_1  v_z_1
    ...                    ...                    ...
    t_1  m_n  R_n  x_n  y_n  z_n  v_x_n  v_y_n  v_z_n
    ...                    ...                    ...
*/
void write_output_data(struct reb_simulation* r) {
    if (do_dry_run) return;

    double t_out = (t_start + r->t * t_unit) / year;
    double m, R, x, y, z, vx, vy, vz;

    FILE* f = fopen(file_out, "a");
    if (!f) error("Failed to open file '%s'", file_out);

    // Header
    if ((is_first_heartbeat) & (!do_restart)) {
        fprintf(f, "# t  m  R  x  y  z  v_x  v_y  v_z  (yr  kg  m  m  m  m/s  m/s  m/s) \n");
    }

    // Append particle data
    for (int i = i_output_start; i < r->N; i++) {
        struct reb_particle p = r->particles[i];

        // Convert to SI
        m = p.m * m_unit;
        R = p.r * l_unit;
        x = p.x * l_unit;
        y = p.y * l_unit;
        z = p.z * l_unit;
        vx = p.vx * v_unit;
        vy = p.vy * v_unit;
        vz = p.vz * v_unit;

        // Output
        fprintf(f, "%e %e %e %e %e %e %e %e %e \n", t_out, m, R, x, y, z, vx, vy, vz);
    }
    fclose(f);
}

/*
    Checkpoint printing and file output every `t_out_step`, and restart files
    every `t_wall_restart`.
*/
void heartbeat(struct reb_simulation* const r) {
    // Change output step
    if ((t_3 > 0) && (t_start + r->t * t_unit > t_3)) {
        t_out_step = t_out_step_3;
    } else if ((t_2 > 0) && (t_start + r->t * t_unit > t_2)) {
        t_out_step = t_out_step_2;
    }

    // Output file and printing
    if (reb_output_check(r, t_out_step / t_unit)) {
        // Revert to the centre-of-mass frame
        if (id_frame_shift >= 0) {
            reb_move_to_com(r);
        }

        reb_output_timing_v2(r, t_end_rel / t_unit);
        write_output_data(r);
    }

    // Activate collisions after some time delay
    if ((collision_mode >= 0) && (r->collision == REB_COLLISION_NONE) &&
        (t_start + r->t * t_unit > t_collision_delay)) {
        r->collision = REB_COLLISION_LINE;
    }

    // Restart file
    if (t_wall_restart > 0.0) {
        if ((!do_dry_run) && (r->walltime > t_wall_restart * restart_step) && (r->t > 0.0)) {
            // Revert to the centre-of-mass frame
            if (id_frame_shift >= 0) {
                reb_move_to_com(r);
            }

            // Save the previous restart file
            if (!is_first_restart) {
                int ret = rename(file_restart, file_restart_prev);
                if (ret != 0)
                    warning(
                        "Failed to rename previous restart file %s to %s", file_restart,
                        file_restart_prev);
            }

            // New restart file
            printf("# Writing restart file... ");
            reb_output_binary(r, file_restart);
            printf("Done \n");
            restart_step++;

            // Mark first restart done
            if (is_first_restart) is_first_restart = 0;
        }
    }

    // Mark first heartbeat done
    if (is_first_heartbeat) is_first_heartbeat = 0;

    // Check wall-clock runtime
    if (r->walltime / 3600.0 > wall_max_h) {
        printf("Wall-clock runtime exceeded, stopping... \n");
        r->status = REB_EXIT_SUCCESS;
    }

    // Shift to the designated particle's (instantaneous) rest frame
    if (reb_output_check(r, t_out_step_1 / t_unit)) {
        if (id_frame_shift >= 0) {
            // Current position and velocity
            double pos_shift[3];
            double vel_shift[3];
            pos_shift[0] = r->particles[id_frame_shift].x;
            pos_shift[1] = r->particles[id_frame_shift].y;
            pos_shift[2] = r->particles[id_frame_shift].z;
            vel_shift[0] = r->particles[id_frame_shift].vx;
            vel_shift[1] = r->particles[id_frame_shift].vy;
            vel_shift[2] = r->particles[id_frame_shift].vz;

            // Shift all particles
            for (int i = 0; i < r->N; i++) {
                r->particles[i].x -= pos_shift[0];
                r->particles[i].y -= pos_shift[1];
                r->particles[i].z -= pos_shift[2];
                r->particles[i].vx -= vel_shift[0];
                r->particles[i].vy -= vel_shift[1];
                r->particles[i].vz -= vel_shift[2];
            }
        }
    }
}

/*
    Update the particle accelerations due to the J2 moment.

    Based on https://rebound.readthedocs.io/en/latest/c_examples/J2/
*/
void force_J2(struct reb_simulation* r) {
    if (oblate_J2 == 0) return;

    const struct reb_particle p_oblate = r->particles[oblate_id];
    const int N = r->N;

#pragma omp parallel for
    for (int i = oblate_id + 1; i < N; i++) {
        const struct reb_particle p = r->particles[i];
        const double sprx = p.x - p_oblate.x;
        const double spry = p.y - p_oblate.y;
        const double sprz = p.z - p_oblate.z;

        const double prx = sprx * cos(-oblate_obliquity) + sprz * sin(-oblate_obliquity);
        const double pry = spry;
        const double prz = -sprx * sin(-oblate_obliquity) + sprz * cos(-oblate_obliquity);
        const double pr2 = prx * prx + pry * pry + prz * prz;

        const double fac =
            3. * r->G * oblate_J2 * p_oblate.m * p_oblate.r * p_oblate.r / 2. / pow(pr2, 3.5);

        const double pax = fac * prx * (prx * prx + pry * pry - 4. * prz * prz);
        const double pay = fac * pry * (prx * prx + pry * pry - 4. * prz * prz);
        const double paz = fac * prz * (3. * (prx * prx + pry * pry) - 2. * prz * prz);

        r->particles[i].ax += pax * cos(oblate_obliquity) + paz * sin(oblate_obliquity);
        r->particles[i].ay += pay;
        r->particles[i].az += -pax * sin(oblate_obliquity) + paz * cos(oblate_obliquity);
    }
}

/*
    Misc temporary testing
*/
void test() {
    printf("## TEST ##\n");

    exit(0);
}

/*
    Print help information.
*/
void print_help() {
    printf(
        "GIHR REBOUND simulations.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "sim : str\n"
        "    The name of the simulation to run. There should exist a folder of this name\n"
        "    with the initial conditions and config files (see load functions) inside it.\n"
        "\n"
        "-h, --help\n"
        "    Print help information and exit.\n"
        "\n"
        "-c <str>, --config=<str>\n"
        "    The name of the config file to load, assumed to be inside the simulation\n"
        "    directory (or a relative path from there). Default `config.txt`.\n"
        "\n"
        "-r, --restart\n"
        "    Restart the simulation from a previously made restart file, `restart.bin`.\n"
        "\n"
        "-n <int>, --n_per_split=<int>\n"
        "-i <int>, --i_split=<int>\n"
        "    For splitting up test particles between runs for crude parallelisation, the\n"
        "    number of test particles per run and the index of this run. e.g. `-n 5 -i 0`\n"
        "    with 25 test particles: simulate only the first 5 test particles, then\n"
        "    `-n 5 -i 1`: simulate the next 5, up to `-n -i 4`: simulate the last 5.\n"
        "\n"
        "-C, --clean\n"
        "    Remove any existing output files for a clean start.\n"
        "\n"
        "-d, --dry\n"
        "    Do a dry run without writing any output files.\n"
        "\n"
        "-s <int>, --steps=<int>\n"
        "    Run for only this many steps.\n"
        "\n"
        "-w <float>, --wall_max_h=<float>\n"
        "    The wall-clock maximum runtime in hours, defaults to 71.5 h.\n"
        "\n"
        "-t, --test\n"
        "    Run temporary testing code.\n");
}

/*
    Parse input parameters (see main docstring).
*/
void parse_parameters(int argc, char** argv) {
    // Short options
    const char* const short_opts = ":hc:rn:i:Cds:w:t";
    // Full options
    const struct option long_opts[] = {
        {"help", no_argument, NULL, 'h'},
        {"config", required_argument, NULL, 'c'},
        {"restart", no_argument, NULL, 'r'},
        {"n_per_split", required_argument, NULL, 'n'},
        {"i_split", required_argument, NULL, 'i'},
        {"clean", no_argument, NULL, 'C'},
        {"dry", no_argument, NULL, 'd'},
        {"steps", required_argument, NULL, 's'},
        {"wall_max_h", required_argument, NULL, 'w'},
        {"test", no_argument, NULL, 't'},
        {0, 0, 0, 0}};

    char file_config_in[999];
    file_config_in[0] = '\0';

    // Parse options
    while (1) {
        int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_help();
                exit(0);
            case 'c':
                sprintf(file_config_in, optarg);
                break;
            case 'r':
                do_restart = 1;
                break;
            case 'n':
                n_per_split = atoi(optarg);
                break;
            case 'i':
                i_split = atoi(optarg);
                break;
            case 'C':
                do_clean = 1;
                break;
            case 'd':
                do_dry_run = 1;
                break;
            case 's':
                steps_end = atoi(optarg);
                break;
            case 'w':
                wall_max_h = atof(optarg);
                break;
            case 't':
                test();
                break;
            case ':':
                error("Option %s requires a value. Run with -h for help.", argv[optind - 1]);
            case '?':
                error("Option %s not recognised. Run with -h for help.", argv[optind - 1]);
        }
    }

    // Other parameters
    sprintf(sim, "%s", argv[optind]);

    // Defaults
    if (file_config_in[0] == '\0') {
        sprintf(file_config_in, "config.txt");
    }

    // Derived parameters
    sprintf(file_config, "%s/%s", sim, file_config_in);
}

/*
    Run the simulation!
*/
int main(int argc, char* argv[]) {
    // Print and parse parameters
    printf("`%s", argv[0]);
    for (int i = 1; i < argc; i++) printf(" %s", argv[i]);
    printf("`\n");
    parse_parameters(argc, argv);

    // Splitting inputs
    if (n_per_split != 0) {
        printf("# Split run [%d] with a subset of %d test particles \n", i_split, n_per_split);
    }

    // Load config parameters
    printf("# Reading config parameters from '%s'\n", file_config);
    int num_picle = load_config_parameters(file_config);

    // Set up a new simulation or a restart
    struct reb_simulation* r;
    if (do_restart) {
        printf("# Loading restart info from '%s'\n", file_restart);

        // Test the file can be read
        FILE* f = fopen(file_restart, "rb");
        if (!f) error("Failed to open file '%s'", file_restart);
        fclose(f);

        // Load the simulation from the restart file
        printf("# Ignore the warning about function pointers, they'll be set before running: ");
        fflush(stdout);
        r = reb_create_simulation_from_binary(file_restart);

        // Back up the previous output files
        char command[999];
        sprintf(command, "cp -fv %s %s.prev", file_out, file_out);
        printf("cp -fv ");
        fflush(stdout);
        system(command);
        sprintf(command, "cp -fv %s %s.prev", file_collisions, file_collisions);
        printf("cp -fv ");
        fflush(stdout);
        system(command);

        // Check if already finished
        if (t_start + r->t * t_unit > 0.999999 * t_end) {
            error(
                "Already reached end time (%.2f yr) before restarting",
                t_start + r->t * t_unit / year);
        }

        // Reset wall-clock time
        r->walltime = 0.0;

        printf(
            "# Run %d particles: %d massive, %d test \n", r->N, r->N_active, r->N - r->N_active);
        printf(
            "# Restart from %.2f to %.2f yr (%.2f yr remaining), output step %.2f yr \n",
            (t_start + r->t * t_unit) / year, t_end / year, (t_end_rel - r->t * t_unit) / year,
            t_out_step / year);
    } else {
        // Load initial conditions
        printf("# Reading initial conditions from '%s'\n", file_init_cond);
        double A1_m[num_picle];
        double A1_R[num_picle];
        double A2_pos[num_picle][3];
        double A2_vel[num_picle][3];
        load_initial_conditions(file_init_cond, num_picle, A1_m, A1_R, A2_pos, A2_vel);
        printf("# Loaded %d particles \n", num_picle);

        // Initialise simulation
        r = reb_create_simulation();

        // Assume all non-test particles are listed first in the initial conditions
        // (or already have zero mass set), override the rest to zero mass
        int i_first_test = num_picle - num_test;

        // Count particles (for this split run)
        int num_picle_run = 0;
        int num_test_run = 0;
        int i_test = 0;

        // Set initial conditions
        for (int i = 0; i < num_picle; i++) {
            struct reb_particle p = {0};

            // Set the particle data (converted to internal units)
            if (i < i_first_test) {
                p.m = A1_m[i] / m_unit;
            } else {
                p.m = 0.0;
            }
            p.r = A1_R[i] / l_unit;
            p.x = A2_pos[i][0] / l_unit;
            p.y = A2_pos[i][1] / l_unit;
            p.z = A2_pos[i][2] / l_unit;
            p.vx = A2_vel[i][0] / v_unit;
            p.vy = A2_vel[i][1] / v_unit;
            p.vz = A2_vel[i][2] / v_unit;

            // Add the particle to the simulation
            if (p.m == 0.0) {
                // Test particle, optionally only add a subset
                if ((n_per_split == 0) || ((i_test >= i_split * n_per_split) &&
                                           (i_test < (i_split + 1) * n_per_split))) {
                    reb_add(r, p);
                    num_test_run++;
                    num_picle_run++;
                }
                i_test++;
            } else {
                // Non-test particle
                reb_add(r, p);
                num_picle_run++;
            }
        }
        if ((n_per_split != 0) && (num_test_run == 0)) {
            error(
                "No test particles selected for n_per_split=%d, i_split=%d", n_per_split, i_split);
        }

        reb_move_to_com(r);

        r->N_active = num_picle_run - num_test_run;
        printf(
            "# Run %d particles: %d massive, %d test \n", num_picle_run, r->N_active,
            num_test_run);
        printf(
            "# Integrate %.2f yr from %.2f to %.2f yr, output step %.2f yr \n", t_end_rel / year,
            t_start / year, t_end / year, t_out_step / year);

        // Set small initial timestep
        r->dt = 1e-5;
    }

    // Variable output step
    if (t_2 > 0) {
        printf("#   After %.2f yr: output step %.2f yr \n", t_2 / year, t_out_step_2 / year);
    }
    if (t_3 > 0) {
        printf("#   After %.2f yr: output step %.2f yr \n", t_3 / year, t_out_step_3 / year);
    }

    // Core functions
    r->heartbeat = heartbeat;
    r->integrator = REB_INTEGRATOR_IAS15;

    // Custom "boundary" option to remove orbiting particles that get too close
    // or too far from the designated primary particle
    if (id_rel_prim >= 0) {
        printf("# Relative distance boundary from particle %d \n", id_rel_prim);
        r->boundary = REB_BOUNDARY_RELATIVE;
        r->id_rel_prim = id_rel_prim;
        r->rel_distance_max = rel_distance_max / l_unit;
    }

    // Oblate potential
    if (oblate_id >= 0) {
        printf("# Oblate potential for particle %d \n", oblate_id);
        r->additional_forces = force_J2;
    }

    // Collisions
    if (collision_mode == -1) {
        r->collision = REB_COLLISION_NONE;
    } else if (collision_mode == 0) {
        r->collision = REB_COLLISION_LINE;
        r->collision_resolve = collision_resolve_write_info;
        printf("# Collisions: write info only \n");
    } else if (collision_mode == 1) {
        r->collision = REB_COLLISION_LINE;
        r->collision_resolve = collision_resolve_hardsphere_and_write;
        r->coefficient_of_restitution = coeff_restitution_fixed;
        printf("# Collisions: hard-sphere (coeff restitution = %.3f) \n", coeff_restitution);
    } else if (collision_mode == 2) {
        r->collision = REB_COLLISION_LINE;
        r->collision_resolve = collision_resolve_merge_and_write;
        printf("# Collisions: merge (note: particle radii not updated) \n");
    } else {
        error("Invalid collision_mode %d", collision_mode);
    }
    if (t_collision_delay > 0) {
        // Disable collisions to start with
        r->collision = REB_COLLISION_NONE;
        printf("# Collisions delayed until %.3f yr \n", t_collision_delay / year);
    }
    if (do_dry_run) {
        printf("# Dry run only, no output \n");
    }
    if (steps_end > 0) {
        printf("# Run only %d steps \n", steps_end);
    }

    // Frame shift
    if (id_frame_shift >= 0) {
        printf("# Recurrent frame shifts to particle %d \n", id_frame_shift);
    }

    // For non-zero split runs, only output test-particle data
    if (i_split != 0) i_output_start = r->N_active;

    double E_initial = reb_tools_energy(r);

    // Go!
    printf("\n");
    if (steps_end > 0) {
        reb_steps(r, steps_end);
    } else {
        reb_integrate(r, t_end_rel / t_unit);
    }
    printf("# Done! (%s)\n", sim);

    // Revert to the centre-of-mass frame
    if (id_frame_shift >= 0) {
        reb_move_to_com(r);
    }

    // Create a final restart file
    if ((t_wall_restart > 0.0) & (!do_dry_run)) {
        // Save the previous restart file
        if (!is_first_restart) {
            int ret = rename(file_restart, file_restart_prev);
            if (ret != 0)
                warning(
                    "Failed to rename previous restart file %s to %s", file_restart,
                    file_restart_prev);
        }

        // New restart file
        printf("# Writing restart file... ");
        reb_output_binary(r, file_restart);
        printf("Done \n");
    }

    double E_final = reb_tools_energy(r);
    printf("# Relative energy change: %.2e. \n", fabs((E_final - E_initial) / E_initial));
    printf(
        "# End time: %.2g. Total steps: %llu. Wall time: ", (t_start + r->t) * t_unit / year,
        r->steps_done);
    if (r->walltime / 60.0 < 3.0) {
        printf("%.3f s. \n", r->walltime);
    } else if (r->walltime / 60.0 < 120.0) {
        printf("%.3g min. \n", r->walltime / 60.0);
    } else {
        printf("%.1f h. \n", r->walltime / 3600.0);
    }

    // Done
    reb_free_simulation(r);
    exit(0);
}
