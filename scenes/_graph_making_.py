from pathlib import Path
import sys
import json
from graph_creation import *


BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
FILENAME = "data.json"

data_zalesak_rotation = [
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic_no_clamped_redistro',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial_no_clamped_redistro',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic.json',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite',
]

data_shear_flow_low = [
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_0_linear',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic_no_clamped_redistro',
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_2_polynomial_no_clamped_redistro',
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_2_polynomial',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE1',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE2',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_0_linear',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_1_cubic',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_2_polynomial',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_3_monotone_hermitn'
]

data_simple_plume_low = [
    '/3_simple_plume/cfl_07/plume_2d_low_conserving_0_linear',
    '/3_simple_plume/cfl_07/plume_2d_low_conserving_1_cubic',
    #'/3_simple_plume/cfl_07/plume_2d_low_conserving_1_cubic_no_clamped_redistro',
    '/3_simple_plume/cfl_07/plume_2d_low_conserving_2_polynomial',
    #'/3_simple_plume/cfl_07/plume_2d_low_conserving_2_polynomial_no_clamped_redistro',
    '/3_simple_plume/cfl_07/plume_2d_low_traditional_RK4_0_linear',
]

data_simple_obstacle_low = [
    '/4_simple_obstacle/cfl_7/obstacle_2d_low_conserving_0_linear',
    #'/4_simple_obstacle/cfl_7/obstacle_2d_low_conserving_1_cubic',
    #'/4_simple_obstacle/cfl_7/obstacle_2d_low_conserving_1_cubic_no_clamped_redistro',
    '/4_simple_obstacle/cfl_7/obstacle_2d_low_conserving_2_polynomial',
    '/4_simple_obstacle/cfl_7/obstacle_2d_low_conserving_2_polynomial_no_clamped_redistro',
    '/4_simple_obstacle/cfl_7/obstacle_2d_low_traditional_RK4_0_linear',
]

data_simple_obstacle_high = [
    '/4_simple_obstacle/cfl_20/obstacle_2d_high_conserving_1_cubic',
    '/4_simple_obstacle/cfl_20/obstacle_2d_high_conserving_1_cubic_local_cfl',
    '/4_simple_obstacle/cfl_20/obstacle_2d_high_conserving_2_polynomial',
    '/4_simple_obstacle/cfl_20/obstacle_2d_high_conserving_2_polynomial_local_cfl',
]

data_plume_2d_high = [
    '/simple_plume_2d_high/simple_plume_2d_high_conserving_0_linear',
    #'/simple_plume_2d_high/simple_plume_2d_high_conserving_0_linear_local_cfl',
    '/simple_plume_2d_high/simple_plume_2d_high_conserving_1_cubic',
    #'/simple_plume_2d_high/simple_plume_2d_high_conserving_1_cubic_local_cfl',
    '/simple_plume_2d_high/simple_plume_2d_high_conserving_2_polynomial',
    #'/simple_plume_2d_high/simple_plume_2d_high_conserving_2_polynomial_local_cfl',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_EE1',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_EE2',
    '/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_0_linear',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_0_linear_local_cfl',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_1_cubic',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_1_cubic_local_cfl',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_2_polynomial',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_2_polynomial_local_cfl',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_3_monotone_hermite',
    #'/simple_plume_2d_high/simple_plume_2d_high_traditional_RK4_3_monotone_hermite_local_cfl',
]

data_plume_3d_high = [
    '/3d_final/simple_plume_3d_high/plume_3d_high_conserving_0_linear',
    '/3d_final/simple_plume_3d_high/plume_3d_high_conserving_0_linear_local_cfl',
    '/3d_final/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_conserving_2_polynomial_local_cfl',
    '/3d_final/simple_plume_3d_high/plume_3d_high_traditional_EE1',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_EE2',
    '/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear',
    '/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_0_linear_local_cfl',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_1_cubic_local_cfl',
    '/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_2_polynomial_local_cfl',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_traditional_RK4_3_monotone_hermite_local_cfl',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_conserving_1_cubic',
    #'/3d_final/simple_plume_3d_high/plume_3d_high_conserving_1_cubic_local_cfl',
]

data = data_simple_obstacle_high
#interested_fields = [["testField", "sum"]] # for fixed vel Field
interested_fields = [["fixed_volume", "sum"]] # for plume

# Resolve full paths
input_paths = []
input_datas = []
input_titles = []
for d in data:
    path = (BASEDIR / d.strip('/') / FILENAME).resolve()
    title = path.parent.name

    input_titles.append(title)
    input_paths.append(path)

    with open(path) as f:
	    input_datas.append(json.load(f))


# Create combined Graph
title = f"combined_{FILENAME.removesuffix(".json")}"
title = "combined_data_nur_vier"
output_path = input_paths[0].parent.parent / f"{title}.pdf"

#General Setup
create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=(12, 6))


#Setup for Shear Flow
""" labels = ["Conserving, linear", "Conserving, polynomial", "EE, linear", "RK4, polynomial"]
output_path = input_paths[0].parent.parent / f"shear_flow_mass_plot.pdf"

create_combined_graph_broken_axis(data_array=input_datas, data_names=labels, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, 
                      figsize=[(12, 3), (12, 1)], y_ranges=[[1960, 2200], [0, 800]], y_gap=0.05, smallTitles=["Total Mass"], yLabels=["Mass"])
 """


#Setup for Zalesak Rotation
""" labels = ["Conserving, linear", "Conserving, polynomial", "EE, linear", "RK4, polynomial"]
output_path = input_paths[0].parent.parent / f"zalesak_rotation_mass_plot.pdf"

create_combined_graph_broken_axis(data_array=input_datas, data_names=labels, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, 
                      figsize=[(12, 3), (12, 1)], y_ranges=[[4306, 4335], [2000, 2900]], y_gap=0.05, smallTitles=["Total Mass"], yLabels=["Mass"])
 """