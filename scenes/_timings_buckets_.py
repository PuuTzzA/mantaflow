import os
import glob
import json
import re
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def extract_index(filename):
    m_all = re.findall(r'(\d+)', filename)
    if not m_all:
        return float('inf')
    return int(m_all[-1])

def gather_timings_from_file(path, categories):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    step = data["step"]
    timings = data["timings"]
    total_ms = data["total_ms"]

    if total_ms is None:
        print("total_ms was None, resorting to summing over entries")
        total_ms = sum(int(x.get("ms", 0)) for x in timings)

    matched = {k: 0 for k in categories}
    matched["other"] = 0
    for entry in timings:
        
        nm = entry.get("name", "")
        ms = int(entry.get("ms", 0))

        found = False
        for k in categories:
            if k in nm:
                matched[k] += ms
                found = True

        if not found and ms < 1e6:
            matched["other"] += ms            

    return step, int(total_ms), matched

def collect_folder_timings(path, categories):

    pattern = os.path.join(path, "timings_*.json")
    files = glob.glob(pattern)

    if not files:
        print(f"No files found in {path}")
        return {}

    files.sort(key=lambda p: extract_index(os.path.basename(p)))

    frame_map = []

    for i, file_path in enumerate(files):
        id, total_ms, matched = gather_timings_from_file(file_path, categories)

        if id != i:
            print("FRAME NOT IN ORDER!")

        pressure_ms = matched[categories[0]]
        advection_ms = matched[categories[1]]
        other_ms = matched["other"]
        #other_ms = max(total_ms - pressure_ms - advection_ms, 0)

        frame_map.append({"pressure_sec": pressure_ms / 1000.0, "advection_sec": advection_ms / 1000.0, "other_sec": other_ms / 1000.0})

    return frame_map

def main():
    base_folder = BASEDIR / BASE_SUBFOLDER

    algorithms_timings = []
    algorithm_names = []

    for (dir, name, advection_algo) in dirs:
        folder = base_folder / dir / "timings"

        categories = [CATEGORY_SOLVE_PRESSURE, advection_algo]
        fmap = collect_folder_timings(str(folder), categories)

        # print(fmap)
        algorithms_timings.append(fmap)
        algorithm_names.append(name)


    algorithm_buckets = {}

    frame_range = FRAME_END - FRAME_START + 1 #plus one because bounds are inclusive
    bucket_amount = int(frame_range / BUCKET_SIZE)

    for index, algorithm_timings in enumerate(algorithms_timings):
        
        advection_total_sum = 0
        pressure_total_sum = 0
        other_total_sum = 0
        steps_amount = 0

        buckets = {}
        for i in range(bucket_amount):

            advection_sum = 0
            pressure_sum = 0
            other_sum = 0
            
            for j in range(BUCKET_SIZE):

                idx = i * BUCKET_SIZE + j

                advection_sum += algorithm_timings[idx]["advection_sec"]
                pressure_sum += algorithm_timings[idx]["pressure_sec"]
                other_sum += algorithm_timings[idx]["other_sec"]

                advection_total_sum += algorithm_timings[idx]["advection_sec"]
                pressure_total_sum += algorithm_timings[idx]["pressure_sec"]
                other_total_sum += algorithm_timings[idx]["other_sec"]
                steps_amount += 1

            buckets[f"{i * BUCKET_SIZE + 1} - {i * BUCKET_SIZE + BUCKET_SIZE}"] = \
                {"advection_average": advection_sum / BUCKET_SIZE, "pressure_average": pressure_sum / BUCKET_SIZE, "other_average": other_sum / BUCKET_SIZE}
        
        print(f"Total averages for {algorithm_names[index]}: (step amount: {steps_amount})")
        print(f"advection: {advection_total_sum / steps_amount}")
        print(f"pressure: {pressure_total_sum / steps_amount}")
        print(f"other: {other_total_sum / steps_amount}")
        print(f"total: {(advection_total_sum + pressure_total_sum + other_total_sum) / steps_amount}")

        algorithm_buckets[algorithm_names[index]] = buckets

    with open(OUTPUT_PATH, "w") as f:
        json.dump(algorithm_buckets, f, indent=4)


# 7_highres_3d_timings
BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()

CATEGORY_SOLVE_PRESSURE = "solvePressure"
CATEGORY_ADVECTION = "simpleSLAdvect"
CATEGORY_ADVECTION_CONSERVING = "massMomentumConservingAdvect"

dirs = [
    ["3d_highres_obstacle_cfl_30_traditional_RK4_0_linear_local_cfl", "not conserving, linear interpolation", CATEGORY_ADVECTION],
    ["3d_highres_obstacle_cfl_30_traditional_RK4_2_polynomial_local_cfl_v2", "not conserving, cubic interpolation", CATEGORY_ADVECTION],
    ["3d_highres_obstacle_cfl_30_conserving_0_linear_local_cfl", "conserving, linear interpolation", CATEGORY_ADVECTION_CONSERVING],
    ["3d_highres_obstacle_cfl_30_conserving_2_polynomial_local_cfl_v3", "conserving, cubic interpolation", CATEGORY_ADVECTION_CONSERVING],
]

BASE_SUBFOLDER = "7_highres_3d_without_vdbs"
OUTPUT_PATH = BASEDIR / BASE_SUBFOLDER / "timings_averages.json"

FRAME_START = 0
FRAME_END = 39
BUCKET_SIZE = 5

if __name__ == "__main__":
    main()