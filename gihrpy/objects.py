"""
Objects for GIHR, accumulated from all objects files.
"""
from gihrpy.objects_phodei import *

# Dictionaries of all of each class's instances so they can be selected by their name
Di_init_prof = InitProf._Di_init_prof
Di_init_cond = InitCond._Di_init_cond
Di_sim = Simulation._Di_simulation
Di_reb_sim = RebSim._Di_reb_sim
Di_sim_set = SimSet._Di_sim_set

# ========
# Main
# ========
if __name__ == "__main__":
    print("\n====  GIHR  objects.py  ====\n")
