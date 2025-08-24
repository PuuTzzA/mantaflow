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

#PATH_TO_MANTA = "/home/tputzer/Documents/BA/mantaflow/manta/"
#PATH_TO_MANTA_VDB = "/home/tputzer/Documents/BA/mantaflow/mantaVDB/"

PATH_TO_MANTA = Path("../../manta").expanduser().resolve()
PATH_TO_MANTA_VDB = Path("../../mantaVDB").expanduser().resolve()

#---------------------------------------------------------------------------------

PATH_TO_SIMPLE_PLUME_SCENE = "../scenes/mass_momentum_conserving.py"
PATH_TO_SIMPLE_PLUME_HIGHRES_SCENE = "../scenes/mass_momentum_conserving_highres.py"
PATH_TO_FIXEL_VELOCITY_SCENE = "../scenes/mmc_fixed_velocity_field.py"
PATH_TO_WATER_SCENE = "../scenes/flip01_mass_momentum_conserving.py"

FIXED_VEL_ZALESAK_ROTATION_PATHS = [
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear.json', False], 
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic.json', False],     
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial_no_clamped_redistro.json', False],           
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial.json', False],           
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1.json', False],
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2.json', False],                   
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear.json', False],          
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic.json', False],         
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial.json', False],      
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite.json', False]
]

FIXED_VEL_SHEAR_FLOW_LOW_PATHS = [
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_conserving_0_linear.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_conserving_1_cubic_no_clamped_redistro.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_conserving_1_cubic.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_conserving_2_polynomial_no_clamped_redistro.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_conserving_2_polynomial.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_EE1.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_EE2.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_RK4_0_linear.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_RK4_1_cubic.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_RK4_2_polynomial.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow_low/shear_flow_traditional_RK4_3_monotone_hermite.json', False]
]

SIMPLE_PLUME_3D_HIGH = [
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_0_linear.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_0_linear_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial.json', True],
    ['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_EE1.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_EE2.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_1_cubic_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_1_cubic.json', True],
]

HIGHRES_PATHS = [
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_conserving_0_linear_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_conserving_2_polynomial_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_traditional_RK4_0_linear_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_traditional_RK4_2_polynomial_local_cfl.json', True],
]

for param, requiresVDB in HIGHRES_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_HIGHRES_SCENE, param])
