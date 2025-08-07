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

PATH_TO_SIMPLE_PLUME_SCENE = "../scenes/mass_momentum_conserving.py"
PATH_TO_FIXEL_VELOCITY_SCENE = "../scenes/mmc_fixed_velocity_field.py"

SIMPLE_PLUME_OBSTACLE_2D_LOW_PATHS = [
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_0_linear.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_1_cubic.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_2_polynomial.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE1.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE2.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_0_linear.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_1_cubic.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_2_polynomial.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_3_monotone_hermite.json', False]
]

SIMPLE_PLUME_OBSTACLE_2D_HIGH_PATHS = [
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_0_linear.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_1_cubic.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_2_polynomial.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE1.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE2.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_0_linear.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_1_cubic.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_2_polynomial.json', False],
    ['../scenes/test_cases/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_3_monotone_hermite.json', False]
]

SIMPLE_PLUME_OBSTACLE_3D_HIGH_PATHS = [
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_0_linear.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_1_cubic.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_2_polynomial.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_EE1.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_EE2.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_0_linear.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_1_cubic.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_2_polynomial.json', True],
    ['../scenes/test_cases/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_3_monotone_hermite.json', True]
]

FIXED_VEL_ZALESAK_ROTATION_PATHS = [
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear.json', False], 
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic.json', False],                
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial.json', False],           
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1.json', False],
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2.json', False],                   
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear.json', False],          
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic.json', False],         
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial.json', False],      
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite.json', False]
]

FIXED_VEL_SHEAR_FLOW_PATHS = [
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_0_linear.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_1_cubic.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_2_polynomial.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_EE1.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_EE2.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_0_linear.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_1_cubic.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_2_polynomial.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_3_monotone_hermite.json', False]
]

for param, requiresVDB in SIMPLE_PLUME_OBSTACLE_3D_HIGH_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_SCENE, param])


""" for param, requiresVDB in FIXED_VEL_ZALESAK_ROTATION_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_FIXEL_VELOCITY_SCENE, param])
 """