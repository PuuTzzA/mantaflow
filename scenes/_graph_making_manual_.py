from pathlib import Path
import sys
import json
from graph_creation import *

FILE_1 = (Path(__file__).parent.parent / "exports/test/oktest/data.json").resolve()
FILE_2 = (Path(__file__).parent.parent / "exports/test/okletsgo/data.json").resolve()

EXPORT_PATH = (Path(__file__).parent.parent / "analysis/data.png").resolve()

data_1 = None
with open(FILE_1) as f:
	data_1 = json.load(f)
data_2 = None
with open(FILE_2) as f:
	data_2 = json.load(f)

create_combined_graph(data_array=[data_1, data_2], data_names=["NOT COSER", "CONSERVING"], interested_fields=[["fixed_volume", "sum"], ["vel_magnitude", "max"]], 
                      title="TEST CREATION", include_cfl=True, include_dt=False, include_extra_stats=False, export_path=EXPORT_PATH)
