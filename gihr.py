"""
Giant Impacts at High Resolution (GIHR) Python main script.

Parameters
----------
func : str
    The name of the function to run. Or instead pass "help" to print the full
    list of available functions and any arguments they need.

--args, -a : str (opt.)
    Any arguments required to be passed to the chosen function. Accepts multiple
    values. Optional arguments can usually be passed as "." to use a default,
    see e.g. check_bool_or_def() in jkeger.py, etc.

--bonus, -b : str (opt.)
    Any bonus options for optional behaviour in the chosen function. Accepts
    multiple values.

--silent, -s
    Redirect standard output to a temporary file, so no output printed to the screen.

--paper, -p
    Tweak plotting bahaviour for making paper figures, e.g. pdf instead of png.

--skip_exist, -k
    Skip making a plot if the file already exists.
"""

from jkeger import *
import gihrpy.utilities as ut
import gihrpy.init_cond as ic
import gihrpy.snapshot_plots as sn
import gihrpy.scenario_plots as sc
import gihrpy.rebound_plots as re
import gihrpy.custom_plots as cu

import argparse

# List of all the functions available to be called by the main program
A1_function = [
    ic.gen_init_prof,
    ic.plot_init_prof,
    ic.gen_init_cond,
    ic.gen_impact_init_cond,
    ic.conv_snap_to_init_cond,
    ic.reb_gen_init_cond,
    #
    sn.plot_snapshots,
    sn.plot_particle_params,
    sn.plot_param_hist_cum,
    #
    sc.plot_sim_param_subsets,
    #
    re.reb_plot_orbit_evol,
    re.reb_plot_pos_trace,
    #
    cu.plot_custom_tile,
    cu.print_results_table,
]
Di_name_function = {func.__name__: func for func in A1_function}


def choose_function(func_choice):
    """Select the function to run from the list of available choices.

    Or, print help text.

    Parameters
    ----------
    func_choice : str
        The name of the function to run.

    Returns
    -------
    func : function
        The function to run.
    """
    # Print help text instead
    if func_choice == "help":
        print("Available functions to run (and all their arguments):")
        # Function names and arguments
        for func in A1_function:
            name = func.__name__
            align = 26
            pad = align + 3

            # Align start of arguments
            string = "  %s(%s)" % (
                add_whitespace(name, align),
                ", ".join(inspect.getfullargspec(func).args),
            )

            # Nice alignment of arguments if longer than one terminal line
            line = 80
            print(string[:line])
            if len(string) > line:
                string = string[line:]
                if string[0] == " ":
                    string = string[1:]
                print(" " * pad + string[: line - pad])
                while len(string) > pad:
                    string = string[line - pad :]
                    print(" " * pad + string[: line - pad])
        print("")
        exit()

    # Return the function with the chosen name
    try:
        return Di_name_function[func_choice]
    except KeyError:
        # If no function matches then print the available list
        print("Unrecognised function choice: ", func_choice)
        print("Available options:")
        for func in A1_function:
            print("  %s()" % func.__name__)
        print(
            "  help                Print the available functions with their arguments."
        )
        exit()


def prep_parser():
    """Prepare the sys args parser."""
    parser = argparse.ArgumentParser()

    # Positional arguments
    parser.add_argument(
        "func", default=None, help="The name of the function to run, or `help`."
    )

    # Optional arguments
    parser.add_argument(
        "-a",
        "--args",
        nargs="*",
        help="The arguments to be passed to the chosen function.",
    )

    parser.add_argument(
        "-b",
        "--bonus",
        nargs="*",
        default=[""],
        help="Any bonus options for optional behaviour in the chosen function.",
    )

    parser.add_argument(
        "-s",
        "--silent",
        action="store_true",
        help="Redirect standard output to a temporary file.",
    )

    parser.add_argument(
        "-p",
        "--paper",
        action="store_true",
        help="Tweak plotting bahaviour, e.g. pdf instead of png.",
    )

    parser.add_argument(
        "-k",
        "--skip_exist",
        action="store_true",
        help="Skip making a plot if the file already exists.",
    )

    return parser


# ========
# Main
# ========
if __name__ == "__main__":
    # Parse arguments
    parser = prep_parser()
    args = parser.parse_args()

    # For globally relevant bool flags, set a variable in utilities.py for easy access
    if args.silent:
        ut.do_silent = True
        # Redirect output
        sys.stdout = open("tmp_silent.txt", "w")
        print("")
    if args.skip_exist:
        ut.do_skip_exist = True
    if args.paper:
        ut.do_paper = True

    # ========
    # Function to run
    # ========
    function = choose_function(args.func)

    # Provided arguments
    A1_arg = args.args

    # Store inputs for global reference
    ut.func_run = args.func
    ut.A1_arg = A1_arg
    ut.A1_bonus = args.bonus

    # ========
    # Run
    # ========
    # In case of any errors, catch to print the inputs after the error message
    try:
        function(*A1_arg)
    except Exception as e:
        print("")
        print(traceback.format_exc())
        msg = "%s(" % function.__name__
        for i_arg, arg in enumerate(A1_arg):
            msg += "%s" % arg
            if i_arg < len(A1_arg) - 1:
                msg += ", "
        msg += ")"
        print("Error for:", msg)

    # Remove temporary matplotlib folder and files
    shutil.rmtree(tmp_dir_mpl)
