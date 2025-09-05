import os
import glob
import json
import re
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def set_custom_margins(ax, left=0, right=0, bottom=0, top=0):
    ax.autoscale(enable=True, axis='both', tight=True)  # no built-in margins

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    dx = x1 - x0
    dy = y1 - y0

    ax.set_xlim(x0 - dx*left,  x1 + dx*right)
    ax.set_ylim(y0 - dy*bottom, y1 + dy*top)

    ax.set_autoscalex_on(False)
    ax.set_autoscaley_on(False)

def main(): 
    with open(INPUT_PATH, 'r', encoding='utf-8') as f:
        data = json.load(f)

    buckets_num = 0
    pressure_arrays = []
    advection_arrays = []
    other_arrays = []

    algorithm_names = []
    bucket_names = []

    first = True
    for algorithm, averages in data.items():
        algorithm_names.append(algorithm)

        a = []
        p = []
        o = []

        for bucket, timings in averages.items():    
            if first:
                buckets_num += 1
                bucket_names.append(bucket)

            for category, sec in timings.items():
                
                if category == "advection_average":
                    a.append(sec)
                elif category == "pressure_average":
                    p.append(sec)
                elif category == "other_average":
                    o.append(sec)

        first = False
        advection_arrays.append(a)
        pressure_arrays.append(p)
        other_arrays.append(o)


    frame_width = BAR_WIDTH * 4 + BAR_PADDING * 3
    fig, ax = plt.subplots(figsize=FIGSIZE, constrained_layout=True)

    x_positions = np.array([i * frame_width + BAR_MARGIN * i for i in range(buckets_num)])
    x_offsets = np.array([BAR_WIDTH + BAR_PADDING for i in range(buckets_num)])

    for k in range(len(algorithm_names)):
        x = x_positions + x_offsets * k

        p = np.array(pressure_arrays[k])
        a = np.array(advection_arrays[k])
        o = np.array(other_arrays[k])

        # solid pressure part
        ax.bar(x, p, width=BAR_WIDTH, color=COLORS[k][0], align='center', linewidth=LINEWIDTH, edgecolor=COLORS[k][0])

        # advection: outline-only (no fill)
        if True or np.any(a > 0):
            ax.bar(x, a, width=BAR_WIDTH, bottom=p, facecolor=COLORS[k][1],
                   edgecolor=COLORS[k][0], linewidth=LINEWIDTH, align='center', hatch='//', alpha=0.05)
            ax.bar(x, a, width=BAR_WIDTH, bottom=p, facecolor="none",
                   edgecolor=COLORS[k][0], linewidth=LINEWIDTH, align='center', hatch='//')


        # other: diagonal hatch (white fill, hatch lines in dir color)
        if True or np.any(o > 0):
            ax.bar(x, o, width=BAR_WIDTH, bottom=p+a, facecolor="none",
                   edgecolor=COLORS[k][0], linewidth=LINEWIDTH, align='center', alpha=1)


            #ax.bar(x, o, width=BAR_WIDTH, bottom=s + a, facecolor='white',
            #       edgecolor=colors[k], linewidth=LINEWIDTH, hatch='O', align='center')

    # x labels are the frame ids
    ax.tick_params(axis='x', length = 0)
    ax.set_xticks(x_positions + BAR_WIDTH * 1.5 + BAR_PADDING * 1.5)
    ax.set_xticklabels(bucket_names)

    ax.set_xlabel("Bucket")
    ax.set_ylabel("Average runtime in seconds")

    # Build legends:
    # 1) colors -> folder mapping
    color_patches = [Patch(facecolor=COLORS[i][0], linewidth=LINEWIDTH, label=algorithm_names[i], edgecolor=COLORS[i][0]) for i in range(len(algorithm_names))]

    # 2) black-style legend to illustrate solid / outline / striped (as requested)
    style_patches = [
        Patch(facecolor='black', edgecolor='black', label=f"pressure projection (solid)"),
        Patch(facecolor='none', edgecolor='black', hatch= '//', linewidth=LINEWIDTH, label=f"advection (striped)"),
        Patch(facecolor='none', edgecolor='black', linewidth=LINEWIDTH , label="other (outline)")
    ]

    # place legends to the right of the plot
    legend1 = ax.legend(handles=color_patches, loc='center left',
                        bbox_to_anchor=BBOX_LEGEND, borderaxespad=0.)
    
    ax.add_artist(legend1)

    legend2 = ax.legend(handles=style_patches, loc='center right',
              bbox_to_anchor=BBBOX_EXPLANATION, borderaxespad=0.)

    #plt.tight_layout(rect=(0, 0, 1, 1))  # leave space on the right for legends

    set_custom_margins(ax, MARGINS[0], MARGINS[1], MARGINS[2], MARGINS[3])
    ax.grid(True, axis="y", linestyle="-", color='gray', linewidth=0.5, alpha=0.3)

    # ensure output dir exists
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    plt.savefig(OUTPUT_PATH, dpi=300,
                bbox_inches='tight',
                bbox_extra_artists=[legend1, legend2],
                pad_inches=0.05)
    print(f"Saved figure to {OUTPUT_PATH} (frames plotted: {buckets_num})")


# 7_highres_3d_timings
BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
BASE_SUBFOLDER = "7_highres_3d_without_vdbs"
INPUT_PATH = BASEDIR / BASE_SUBFOLDER / "timings_averages.json"
OUTPUT_PATH = BASEDIR / BASE_SUBFOLDER / "timings.pdf"

FIGSIZE = (15, 6.5)
LINEWIDTH = 1
BAR_WIDTH = 0.8
BAR_PADDING = 0.2
BAR_MARGIN = 1

offset = 0.0585
BBOX_LEGEND = (offset, 1)
BBBOX_EXPLANATION = (1 - offset, 1)
MARGINS = (0.03, 0.03, 0, 0.3) #left, right, bottom, top

blue = plt.cm.tab20.colors[0]
light_blue = blue
orange = plt.cm.tab20.colors[2]
light_orange = orange
green = plt.cm.tab20.colors[4]
light_green = green
violet = plt.cm.tab20.colors[8]
light_violet = violet

#blue = plt.cm.tab20.colors[0]
#light_blue = plt.cm.tab20c.colors[2]
#orange = plt.cm.tab20.colors[2]
#light_orange = plt.cm.tab20c.colors[7]
#green = plt.cm.tab20.colors[4]
#light_green = plt.cm.tab20c.colors[10]
#violet = plt.cm.tab20.colors[8]
#light_violet = plt.cm.tab20c.colors[14]
#light_blue = light_orange = light_green = light_violet = "none"

COLORS = [[violet, light_violet], [green, light_green], [blue, light_blue], [orange, light_orange]]

allTextSize = 24

plt.style.use("classic")

plt.rcParams["legend.fontsize"] = allTextSize  # or a number like 10
plt.rcParams["axes.labelsize"] = allTextSize
plt.rcParams["xtick.labelsize"] = allTextSize - 2
plt.rcParams["ytick.labelsize"] = allTextSize - 2
mpl.rcParams['hatch.linewidth'] = 1  # previous pdf hatch linewidth

mpl.rcParams.update({
    "text.usetex": False,
    "font.family": "TeX Gyre Pagella",
    "mathtext.fontset": "cm",   # uses Matplotlib’s Computer Modern lookalike
    "axes.xmargin": 0,
    })

if __name__ == "__main__":
    main()