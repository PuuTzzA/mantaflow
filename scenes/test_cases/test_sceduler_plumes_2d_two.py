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
PATH_TO_FIXEL_VELOCITY_SCENE = "../scenes/mmc_fixed_velocity_field.py"
PATH_TO_WATER_SCENE = "../scenes/flip01_mass_momentum_conserving.py"
PATH_TO_SIMPLE_PLUME_HIGHRES_SCENE = "../scenes/mass_momentum_conserving_highres.py"

FIXED_VEL_ZALESAK_ROTATION_PATHS = [
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear.json', False], 
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic.json', False],     
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial_no_clamped_redistro.json', False],           
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial.json', False],           
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1.json', False],
    #['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2.json', False],                   
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear.json', False],          
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic.json', False],         
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial.json', False],      
    ['../scenes/test_cases/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite.json', False]
]

FIXED_VEL_SHEAR_FLOW_PATHS = [
#    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_0_linear.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_1_cubic.json', False],
#    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_2_polynomial_no_clamped_redistro.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_conserving_2_polynomial.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_EE1.json', False],
    #['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_EE2.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_0_linear.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_1_cubic.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_2_polynomial.json', False],
    ['../scenes/test_cases/fixed_vel_shear_flow/shear_flow_traditional_RK4_3_monotone_hermite.json', False]
]

PLUME_2D_LOW_PATHS = [
#    ['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_0_linear.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_1_cubic.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_1_cubic_local_cfl_no_clamped_redistro.json', False],
#    ['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_2_polynomial.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_2_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_2_polynomial_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_conserving_2_polynomial_local_cfl_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_EE1.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_EE2.json', False],
#    ['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_0_linear.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_1_cubic.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_1_cubic_local_cfl.json', False],
    ['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_2_polynomial.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_3_monotone_hermite.json', False],
    #['../scenes/test_cases/simple_plume_2d_low/plume_2d_low_traditional_RK4_3_monotone_hermite_local_cfl.json', False],
]

PLUME_2D_HIGH_PATHS = [
    ['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_0_linear.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_1_cubic.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_1_cubic_local_cfl_no_clamped_redistro.json', False],
#    ['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_2_polynomial.json', False],
#    ['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_2_polynomial_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_conserving_2_polynomial_local_cfl_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_EE1.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_EE2.json', False],
    ['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_0_linear.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_1_cubic.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_2_polynomial.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_3_monotone_hermite.json', False],
    #['../scenes/test_cases/simple_plume_2d_high/plume_2d_high_traditional_RK4_3_monotone_hermite_local_cfl.json', False],
]

OBSTACLE_2D_LOW_PATHS = [
#    ['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_0_linear.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_1_cubic.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_1_cubic_local_cfl_no_clamped_redistro.json', False],
#    ['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_2_polynomial.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_2_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_2_polynomial_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_conserving_2_polynomial_local_cfl_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_EE1.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_EE2.json', False],
#    ['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_0_linear.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_1_cubic.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_1_cubic_local_cfl.json', False],
    ['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_2_polynomial.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_3_monotone_hermite.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_low/obstacle_2d_low_traditional_RK4_3_monotone_hermite_local_cfl.json', False],
]

OBSTACLE_2D_HIGH_PATHS = [
#    ['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_0_linear.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_1_cubic.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_1_cubic_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_1_cubic_local_cfl_no_clamped_redistro.json', False],
#    ['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_2_polynomial.json', False],
    ['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_2_polynomial_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_conserving_2_polynomial_local_cfl_no_clamped_redistro.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_EE1.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_EE2.json', False],
    ['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_0_linear.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_0_linear_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_1_cubic.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_1_cubic_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_2_polynomial.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_2_polynomial_local_cfl.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_3_monotone_hermite.json', False],
    #['../scenes/test_cases/simple_obstacle_2d_high/obstacle_2d_high_traditional_RK4_3_monotone_hermite_local_cfl.json', False],
]

HIGHRES_PATHS = [
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_40_conserving_2_polynomial_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_40_conserving_0_linear_local_cfl.json', True],

    #['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_20_conserving_0_linear_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_20_conserving_2_polynomial_local_cfl.json', True],

    #['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_10_conserving_0_linear_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_10_conserving_2_polynomial_local_cfl.json', True],

    ['../scenes/test_cases/3d_highres/3d_highres_no_obstacle_cfl_40_conserving_2_polynomial_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_no_obstacle_cfl_40_conserving_0_linear_local_cfl.json', True],
]

SIMPLE_PLUME_3D_HIGH = [
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_0_linear.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_0_linear_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_1_cubic.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_1_cubic_local_cfl.json', True],
#    ['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial_local_cfl.json', True],
    ['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_EE1.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_EE2.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic_local_cfl.json', True],
#    ['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial.json', True],
#    ['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial_local_cfl.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite.json', True],
    #['../scenes/test_cases/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite_local_cfl.json', True],
]

DIFFERENT_CFL_2D_PATHS = [
#    # Same Resolution 200
#    ['../scenes/test_cases/different_cfl_2d_res_200/40_dif_cfl_conserving_polynomial_local_cfl.json', False],                        
#    ['../scenes/test_cases/different_cfl_2d_res_200/30_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_200/20_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_200/10_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_200/05_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_200/03_dif_cfl_conserving_polynomial_local_cfl.json', False],

#    # Same Resolution 200
#    ['../scenes/test_cases/different_cfl_2d_res_150/40_dif_cfl_conserving_polynomial_local_cfl.json', False],                        
#    ['../scenes/test_cases/different_cfl_2d_res_150/30_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_150/20_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_150/10_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_150/05_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d_res_150/03_dif_cfl_conserving_polynomial_local_cfl.json', False],


    # Different cfl and resolution
#    ['../scenes/test_cases/different_cfl_2d/02_dif_cfl_conserving_polynomial_local_cfl.json', False],
    ['../scenes/test_cases/different_cfl_2d/02_dif_cfl_traditional_polynomial_local_cfl.json', False],

    ['../scenes/test_cases/different_cfl_2d/05_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d/05_dif_cfl_traditional_polynomial_local_cfl.json', False],

#    ['../scenes/test_cases/different_cfl_2d/40_dif_cfl_conserving_polynomial_local_cfl.json', False],
    ['../scenes/test_cases/different_cfl_2d/40_dif_cfl_traditional_polynomial_local_cfl.json', False],
                          
    ['../scenes/test_cases/different_cfl_2d/30_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d/30_dif_cfl_traditional_polynomial_local_cfl.json', False],
                           
#    ['../scenes/test_cases/different_cfl_2d/20_dif_cfl_conserving_polynomial_local_cfl.json', False],
    ['../scenes/test_cases/different_cfl_2d/20_dif_cfl_traditional_polynomial_local_cfl.json', False],
                          
    ['../scenes/test_cases/different_cfl_2d/60_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d/60_dif_cfl_traditional_polynomial_local_cfl.json', False],

#    ['../scenes/test_cases/different_cfl_2d/80_dif_cfl_conserving_polynomial_local_cfl.json', False],
    ['../scenes/test_cases/different_cfl_2d/80_dif_cfl_traditional_polynomial_local_cfl.json', False],

    ['../scenes/test_cases/different_cfl_2d/10_dif_cfl_conserving_polynomial_local_cfl.json', False],
#    ['../scenes/test_cases/different_cfl_2d/10_dif_cfl_traditional_polynomial_local_cfl.json', False],
]

HIGHRES_PATHS = [
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_conserving_2_polynomial_local_cfl.json', True],
    # ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_traditional_RK4_2_polynomial_local_cfl.json', True],

    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_conserving_0_linear_local_cfl.json', True],
    ['../scenes/test_cases/3d_highres/3d_highres_obstacle_cfl_30_traditional_RK4_0_linear_local_cfl.json', True],
]

WATER_PATHS = [
    #['../scenes/test_cases/simple_water/simple_water_dam_FLIP_traditional.json', False],
    #['../scenes/test_cases/simple_water/simple_water_dam_NO_FLIP_conserving.json', False],
    #['../scenes/test_cases/simple_water/simple_water_dam_NO_FLIP_traditional.json', False],
    #['../scenes/test_cases/simple_water/simple_water_drop_FLIP_traditional.json', False],
    ['../scenes/test_cases/simple_water/simple_water_drop_NO_FLIP_conserving.json', False],
    ['../scenes/test_cases/simple_water/simple_water_drop_NO_FLIP_traditional.json', False],
]

""" 
for param, requiresVDB in FIXED_VEL_SHEAR_FLOW_PATHS + FIXED_VEL_ZALESAK_ROTATION_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_FIXEL_VELOCITY_SCENE, param])

for param, requiresVDB in PLUME_2D_LOW_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_SCENE, param])

os.chdir(PATH_TO_MANTA)
subprocess.run(["make", "-j4"]) 

for param, requiresVDB in OBSTACLE_2D_LOW_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_SCENE, param])

 """

for param, requiresVDB in DIFFERENT_CFL_2D_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_SCENE, param])


""" for param, requiresVDB in SIMPLE_PLUME_3D_HIGH:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_SCENE, param])
 """

""" for param, requiresVDB in HIGHRES_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_SIMPLE_PLUME_HIGHRES_SCENE, param])
 """


""" for param, requiresVDB in WATER_PATHS:

    if requiresVDB:
        os.chdir(PATH_TO_MANTA_VDB)
    else:
        os.chdir(PATH_TO_MANTA)

    subprocess.run(["./manta", PATH_TO_WATER_SCENE, param]) 
 """