import os
import glob
import json
import re
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()   

def extract_index(filename):
    """Return integer index found in filename (last group of digits), or inf if none."""
    m_all = re.findall(r'(\d+)', filename)
    if not m_all:
        return float('inf')
    return int(m_all[-1])


def gather_timings_from_file(path, name_keys):
    """
    Read one JSON file and return (step, total_ms, dict of matched_ms for each key).
    """
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    step = data.get("step", None)
    timings = data.get("timings", [])
    total_ms = data.get("total_ms")
    if total_ms is None:
        total_ms = sum(int(x.get("ms", 0)) for x in timings)

    matched = {k: 0 for k in name_keys}
    for entry in timings:
        nm = entry.get("name", "")
        ms = int(entry.get("ms", 0))
        for k in name_keys:
            if k in nm:
                matched[k] += ms

    return step, int(total_ms), matched


def main():
    pattern = os.path.join(FOLDER, "timings_*.json")
    files = glob.glob(pattern)
    if not files:
        print(f"No files found in {FOLDER}")
        return

    files.sort(key=lambda p: extract_index(os.path.basename(p)))

    name_keys = [CATEGORY_SOLVE_PRESSURE, CATEGORY_ADVECTION]

    steps = []
    solve_list, advect_list, other_list = [], [], []
    total_ms_list = []

    for i, fpath in enumerate(files):
        step, total_ms, matched = gather_timings_from_file(fpath, name_keys)

        step_label = step if step is not None else i
        steps.append(step_label)

        solve_ms = matched[CATEGORY_SOLVE_PRESSURE]
        advect_ms = matched[CATEGORY_ADVECTION]
        other_ms = total_ms - solve_ms - advect_ms
        if other_ms < 0:
            other_ms = 0

        solve_list.append(solve_ms)
        advect_list.append(advect_ms)
        other_list.append(other_ms)
        total_ms_list.append(total_ms)

    n = len(solve_list)

    solve_list = np.array(solve_list  , dtype=np.float64)
    advect_list = np.array(advect_list, dtype=np.float64)
    other_list = np.array(other_list  , dtype=np.float64)

    solve_list /= (1000)
    advect_list /= (1000)
    other_list /= (1000)

    avg_solve = np.mean(solve_list)
    avg_advect = np.mean(advect_list)
    avg_other = np.mean(other_list)

    print(f"Frames parsed: {n}")
    print(f"Average solvePressure : {avg_solve:.3f} s")
    print(f"Average advect        : {avg_advect:.3f} s")
    print(f"Average other         : {avg_other:.3f} s")
    print(f"Average total         : {(avg_solve + avg_advect + avg_other):.3f} s")

    # plotting
    x = list(range(n))

    fig, ax = plt.subplots(figsize=(max(10, n * 0.15), 6))
    ax.bar(x, solve_list, label="solvePressure")
    ax.bar(x, advect_list, bottom=solve_list, label="massMomentumConservingAdvect")
    bottom = [s + a for s, a in zip(solve_list, advect_list)]
    ax.bar(x, other_list, bottom=bottom, label="other")

    ax.set_xlabel("Frame")
    ax.set_ylabel("Duration (ms)")
    ax.set_title("Per-frame timings (stacked)")

    # averages
    ax.axhline(avg_solve, linestyle="--", linewidth=1,
               label=f"avg solve: {avg_solve:.2f} s")
    ax.axhline(avg_solve + avg_advect, linestyle="-.", linewidth=1,
               label=f"avg advect: {avg_advect:.2f} s")
    ax.axhline(avg_solve + avg_advect + avg_other, color="gray", linestyle=":",
               linewidth=1, label=f"avg total: {(avg_solve+avg_advect+avg_other):.2f} s")

    ax.legend(loc="upper right")
    plt.tight_layout()

    if OUTFILE:
        plt.savefig(OUTFILE, dpi=150)
        print(f"Saved figure to {OUTFILE}")
    else:
        plt.show()


FOLDER = BASEDIR /  "7_highres_3d/3d_highres_obstacle_cfl_30_traditional_RK4_2_polynomial_local_cfl/timings"
OUTFILE = BASEDIR / "7_highres_3d/3d_highres_obstacle_cfl_30_traditional_RK4_2_polynomial_local_cfl/timings_plot.png"
# ----------------

CATEGORY_SOLVE_PRESSURE = "solvePressure"
CATEGORY_ADVECTION = "simpleSLAdvect"
#CATEGORY_ADVECTION = "massMomentumConservingAdvect"

if __name__ == "__main__":
    main()
