import json
from pathlib import Path
from itertools import product
import shutil

# Base directory to store JSONs
CONTAINER_DIR = "simple_plume_3d_high_highres"
BASE_DIR = (Path(__file__).parent / CONTAINER_DIR).resolve()
shutil.rmtree(BASE_DIR, ignore_errors=True)
BASE_DIR.mkdir(parents=True, exist_ok=True)

# Constant settings
BASE_TITLE = "plume_3d_high_highres"
dimension = 3
resolutionX = 256
resolutionY = 512
resolutionZ = 256

doOpen = True
doObstacle = True

exportData = True
exportImages = False
exportVideos = False
exportVDBs = True  # This is the flag we will output as True/False in the printed list

max_time = 140
maxCFL = 50
dt = 6

# Interpolation method names
interpolation_method_names = {
    -1: "",
    0: "0_linear",
    1: "1_cubic",
    2: "2_polynomial",
    3: "3_monotone_hermite",
}

tracing_function_names = {
    0 : "",
    1 : "local_cfl",
}

# Settings to vary
conserving_options = [True, False]
tracing_methods = ["EE1", "EE2", "RK4"]
interpolation_methods = [0, 1, 2, 3]  # 3 is invalid when conserving=True
tracing_functions = [0, 1]

generated_paths = []

for doConserving in conserving_options:
    if doConserving:
        valid_interps = [0, 1, 2]  # monotoneHermite not allowed
        for interp in valid_interps:
            for tracing_function in tracing_functions:
                interp_name = interpolation_method_names[interp]
                tracing_function_name = tracing_function_names[tracing_function]
                title = f"{BASE_TITLE}_conserving{"_" if interp_name != "" else ""}{interp_name}{"_" if tracing_function != 0 else ""}{tracing_function_name}"
                config = {
                    "title": title,
                    "dimension": dimension,
                    "resolutionX": resolutionX,
                    "resolutionY": resolutionY,
                    "resolutionZ": resolutionZ,
                    "doOpen": doOpen,
                    "doObstacle": doObstacle,
                    "doConserving": True,
                    "tracingMethod": "RK4",  # ignored in this case
                    "tracingFunction": tracing_function,
                    "interpolationMethod": interp,
                    "exportData": exportData,
                    "exportImages": exportImages,
                    "exportVideos": exportVideos,
                    "exportVDBs": exportVDBs,
                    "max_time": max_time,
                    "maxCFL": maxCFL,
                    "dt": dt,
                }

                path = BASE_DIR / f"{title}.json"
                with open(path, "w") as f:
                    json.dump(config, f, indent=2)

                # Store relative path with exportVDBs flag
                relative_path = f"../scenes/test_cases/{CONTAINER_DIR}/{title}.json"
                generated_paths.append([relative_path, exportVDBs])

    else:
        for tracing in tracing_methods:
            if tracing == "RK4":
                interps = interpolation_methods  # all interpolation methods
            else:
                interps = [-1] 
            for interp in interps:
                for tracing_function in tracing_functions:
                    interp_name = interpolation_method_names[interp]
                    tracing_function_name = tracing_function_names[tracing_function]
                    title = f"{BASE_TITLE}_traditional_{tracing}{"_" if interp_name != "" else ""}{interp_name}{"_" if tracing_function != 0 else ""}{tracing_function_name}"
                    config = {
                        "title": title,
                        "dimension": dimension,
                        "resolutionX": resolutionX,
                        "resolutionY": resolutionY,
                        "resolutionZ": resolutionZ,
                        "doOpen": doOpen,
                        "doObstacle": doObstacle,
                        "doConserving": False,
                        "tracingMethod": tracing,
                        "interpolationMethod": interp,
                        "tracingFunction": tracing_function,
                        "exportData": exportData,
                        "exportImages": exportImages,
                        "exportVideos": exportVideos,
                        "exportVDBs": exportVDBs,
                        "max_time": max_time,
                        "maxCFL": maxCFL,
                        "dt": dt,
                    }

                    path = BASE_DIR / f"{title}.json"
                    with open(path, "w") as f:
                        json.dump(config, f, indent=2)

                    # Store relative path with exportVDBs flag
                    relative_path = f"../scenes/test_cases/{CONTAINER_DIR}/{title}.json"
                    generated_paths.append([relative_path, exportVDBs])

# Output the list
for entry in generated_paths:
    print(entry)


