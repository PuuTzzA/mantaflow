# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# You have to run this file from the ubuntu terminal, bc vscode was installed with 
# snap and manta uses DLLs that are not compatible with Snaps isolated environment
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import os
from pathlib import Path
import subprocess
import shutil
import json
import numpy as np
import matplotlib.pyplot as plt
import math

PATH_TO_MANTA = "/home/tputzer/Documents/BA/mantaflow/manta/"
PATH_TO_MANTA_VDB = "/home/tputzer/Documents/BA/mantaflow/mantaVDB/"

#---------------------------------------------------------------------------------

PATH_TO_SCENE = "../scenes/mass_momentum_conserving.py"
PATHS_TO_PARAMS_GAS_3d = [
    #["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_sl_cfl_1.json", True], 
    #["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_sl_cfl_15.json", True], 
    #["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_sl_cfl_30.json", True], 
    ["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_conserving_cfl_1.json", True], 
    ["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_conserving_cfl_15.json", True], 
    ["../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_conserving_cfl_30.json", True], 

    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_sl_cfl_1.json", True], 
    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_sl_cfl_15.json", True], 
    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_sl_cfl_30.json", True], 
    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_conserving_cfl_1.json", True], 
    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_conserving_cfl_15.json", True], 
    ["../scenes/test_cases/simple_obstacle_tests/params_simple_obstacle_3d_conserving_cfl_30.json", True], 
]
                   

for param, requiresVDB in PATHS_TO_PARAMS_GAS_3d:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SCENE, param])