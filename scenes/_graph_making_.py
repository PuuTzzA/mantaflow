from pathlib import Path
import sys
import json
from graph_creation import *


BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
FILENAME = "data.json"

videos_zalesak_rotation = [
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite/',
]

videos = videos_zalesak_rotation

# Resolve full paths
input_paths = []
input_datas = []
input_titles = []
for v in videos:
    path = (BASEDIR / v.strip('/') / FILENAME).resolve()
    title = path.parent.name

    input_titles.append(title)
    input_paths.append(path)

    with open(path) as f:
	    input_datas.append(json.load(f))


# Create combined Graph
title = f"combined_{FILENAME.removesuffix(".json")}"
output_path = input_paths[0].parent.parent / f"combined_{title}.png"

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=[["testField", "sum"]], 
                      title=title, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=(12, 12))
