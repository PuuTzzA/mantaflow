from pathlib import Path
import sys
import json
from graph_creation import *


BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
FILENAME = "data.json"

data_zalesak_rotation = [
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear/',
    #'/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial/',
    #'/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1/',
    #'/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear/',
    #'/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial/',
    #'/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite/',
]

data_shear_flow = [
    '/fixed_vel_shear_flow/shear_flow_conserving_0_linear',
    #'/fixed_vel_shear_flow/shear_flow_conserving_1_cubic',
    '/fixed_vel_shear_flow/shear_flow_conserving_2_polynomial',
    #'/fixed_vel_shear_flow/shear_flow_traditional_EE1',
    #'/fixed_vel_shear_flow/shear_flow_traditional_EE2',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_0_linear',
    #'/fixed_vel_shear_flow/shear_flow_traditional_RK4_1_cubic',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_2_polynomial',
    #'/fixed_vel_shear_flow/shear_flow_traditional_RK4_3_monotone_hermite',
]

data_plume_2d_low = [
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_0_linear',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_1_cubic',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_2_polynomial',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE1',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE2',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_0_linear',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_1_cubic',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_2_polynomial',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_3_monotone_hermite',
]

data_plume_2d_high = [
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_0_linear',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_1_cubic',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_2_polynomial',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE1',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE2',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_0_linear',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_1_cubic',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_2_polynomial',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_3_monotone_hermite',
]

data_plume_3d_high = [
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_0_linear',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_1_cubic',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_conserving_2_polynomial',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_EE1',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_EE2',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_0_linear',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_1_cubic',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_2_polynomial',
    '/simple_plume_obstacle_3d_high/simple_plume_obstacle_3d_high_traditional_RK4_3_monotone_hermite',
]

data = data_shear_flow
interested_fields = [["testField", "sum"]] # for fixed vel Field
#interested_fields = [["fixed_volume", "sum"]] # for plume

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
title = "combined_data_3"
output_path = input_paths[0].parent.parent / f"{title}.pdf"

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=(12, 6))
