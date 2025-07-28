import os
from pathlib import Path
import subprocess
import shutil
import json
import numpy as np
import matplotlib.pyplot as plt
import math


PATH_TO_TEST_CASES = "/home/tputzer/Documents/BA/mantaflow/scenes/test_cases/3d_tests_simple_plume/"
PATH_TO_MANTA = "/home/tputzer/Documents/BA/mantaflow/manta/"
PATH_TO_MANTA_VDB = "/home/tputzer/Documents/BA/mantaflow/mantaVDB/"

subprocess.run(["ls", "ls"])