from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import math

COLOR_THEMES = [
    {
        "main": "blue",
        "mean": "purple",
        "median": "deeppink",
        "min": "darkslategray",
        "max": "darkslategray"
    },
    {
        "main": "red",
        "mean": "darkgoldenrod",
        "median": "sienna",
        "min": "black",
        "max": "black"
    },
    {
        "main": "green",
        "mean": "seagreen",
        "median": "mediumseagreen",
        "min": "darkgreen",
        "max": "darkgreen"
    }
]

def create_combined_graph(data_array, data_names, interested_fields, title, include_cfl=True, include_extra_stats=True, export_path="./exports/combinded.png"):
    frames_sets = []

    for data in data_array:
        frames = data["frame_data"]
        if not frames:
            print("No frame data collected")
            return
        
        field_frames = {}
        if include_cfl:
            field_frames["cfl"] = []
        for interested_field in interested_fields:
            field_frames[f"{interested_field} Sum"] = []

        for key, grid in frames.items():
            if include_cfl:
                field_frames["cfl"].append(grid["cfl"])

            for interested_field in interested_fields:
                field_frames[f"{interested_field} Sum"].append(grid[interested_field]["sum"])
        
        frames_sets.append(field_frames)

    amount = len(interested_fields) + (1 if include_cfl else 0)
    fig, ax = plt.subplots(amount, 1, figsize=(12, 4 * amount), sharex=True)
    ax = np.atleast_1d(ax)

    current_ax = 0
    for key, value in frames_sets[0].items():

        for i in range(len(data_array)):     

            frames = np.array(frames_sets[i][key])
            minimum = np.min(frames)
            maximum = np.max(frames)
            mean = np.mean(frames)
            std_dev = np.std(frames)
            median = np.median(frames)

            data_name = f", {data_names[i]}" if data_names[i] != "" else ""
            ax[current_ax].plot(frames, linestyle='-', linewidth=2.5, color=COLOR_THEMES[i]["main"], label=f'{key}{data_name}', zorder=2)
            if include_extra_stats:
                ax[current_ax].axhline(mean, color=COLOR_THEMES[i]["mean"], linestyle='--', linewidth=1.2, label=f'Mean{data_name}: {mean:.2f}', zorder=1)
                ax[current_ax].axhline(median, color=COLOR_THEMES[i]["median"], linestyle='--', linewidth=1.2, label=f'Median{data_name}: {median:.2f}', zorder=1)
                ax[current_ax].axhline(minimum, color=COLOR_THEMES[i]["min"], linestyle=':', linewidth=1.2, label=f'Min{data_name}: {minimum:.2f}', zorder=1)
                ax[current_ax].axhline(maximum, color=COLOR_THEMES[i]["max"], linestyle=':', linewidth=1.2, label=f'Max{data_name}: {maximum:.2f}', zorder=1)
            ax[current_ax].set_ylabel(key)
            ax[current_ax].set_title(f"{key} Over Time")
            ax[current_ax].legend(loc='best')
            ax[current_ax].grid(True)
            ax[current_ax].ticklabel_format(style='plain', axis='y', useOffset=False)
            
        current_ax += 1

    plt.xlabel("Timesteps")
    fig.suptitle(title, fontsize="x-large")
    plt.tight_layout()
    plt.savefig(export_path, dpi=300)

#NAME = "shear_flow"
#data_1 = None
#data_2 = None
#with open(f"./exports/experiments/{NAME}_Traditional_stats/data.json") as f:
#	data_1 = json.load(f)
#with open(f"./exports/experiments/{NAME}_Conserving_stats/data.json") as f:
#	data_2 = json.load(f)
#
#create_combined_graph(data_array=[data_1, data_2], data_names=["Traditional", "Mass Momemtum Conserving"], interested_fields=["testField"], 
#                      title="Fixed Shear Flow Field", include_cfl=True, include_extra_stats=False)