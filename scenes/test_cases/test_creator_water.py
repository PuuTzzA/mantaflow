import json
from pathlib import Path
from itertools import product
import shutil

# Base directory to store JSONs
CONTAINER_DIR = "simple_water"
BASE_DIR = (Path(__file__).parent / CONTAINER_DIR).resolve()
shutil.rmtree(BASE_DIR, ignore_errors=True)
BASE_DIR.mkdir(parents=True, exist_ok=True)

# Constant settings
BASE_TITLE = "simple_water"
dimension = 2
resolutionX = 64
resolutionY = 64
resolutionZ = 64

layout_options = ["dam", "drop"]
doPLS = True

exportData = True
exportImages = True
exportVideos = True
exportVDBs = False  # This is the flag we will output as True/False in the printed list

max_time = 150
maxCFL = 0.5
dt = 1

# Interpolation method names
interpolation_method_names = {
    -1: "",
    0: "0_linear",
    1: "1_cubic",
    2: "2_polynomial",
}

tracing_function_names = {
    0 : "",
    1 : "local_cfl",
}

flip_not_flip = [True, False]

# Settings to vary
conserving_options = [True, False]
interpolation_methods = [0, 1, 2]  # 3 is invalid when conserving=True
#tracing_functions = [0, 1]

generated_paths = []

for layout in layout_options:
    for doFlip in flip_not_flip:
        conserving_options = [False] if doFlip else [True, False]
        
        for doConserving in conserving_options: 
            #    tracing_function_name = tracing_function_names[tracing_function]
            title = f"{BASE_TITLE}_{layout}_{"FLIP" if doFlip else "NO_FLIP"}_{"conserving" if doConserving else "traditional"}"
            config = {
                "title": title,
                "dimension": dimension,
                "resolutionX": resolutionX,
                "resolutionY": resolutionY,
                "resolutionZ": resolutionZ,
                "layout" : layout,
                "doParticleLevelSet" : doPLS,
                "doFLIP": doFlip,
                "doConserving": doConserving,
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


