from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
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


def create_combined_graph(data_array, data_names, interested_fields, title, include_title=True, include_cfl=True, include_dt=True, include_extra_stats=True, export_path="./exports/combinded.png", figsize=(12,4)):
    """
    Creates and saves a combined multi-line graph showing the evolution of various fields over time.

    Parameters:
        data_array (list of dict): List of datasets, each containing "frame_data" with field values per timestep.
        data_names (list of str): Names corresponding to each dataset for labeling in the graph legend.
        interested_fields (list of tuples): List of (field_name, target_key (like sum, max, ...)) pairs to extract from each frame.
        title (str): The title for the entire figure.
        include_cfl (bool): If True, includes the CFL field in the plots.
        include_dt (bool): If True, includes the timestep (dt) in the plots.
        include_extra_stats (bool): If True, overlays mean, median, min, and max lines for each dataset.
        export_path (str): File path where the final plot image will be saved.
        figsize (tuple of number): Size of one individual graph.
    """
    frames_sets = []

    for data in data_array:
        frames = data["frame_data"]
        if not frames:
            print("No frame data collected")
            return
        
        field_frames = {}
        if include_cfl:
            field_frames["cfl"] = []
        if include_dt:
            field_frames["dt"] = []
        for interested_field, target in interested_fields:
            field_frames[f"{interested_field} {target}"] = []
              
        for key, grid in frames.items():
            if include_cfl:
                field_frames["cfl"].append(grid["cfl"])
            
            if include_dt:
                field_frames["dt"].append(grid["dt"])

            for interested_field, target in interested_fields:
                field_frames[f"{interested_field} {target}"].append(grid[interested_field][target])
        
        frames_sets.append(field_frames)

    #plt.style.use("seaborn-v0_8-deep")
    plt.style.use("seaborn-v0_8-notebook")
    #plt.style.use("default")

    mpl.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "mathtext.fontset": "cm"   # uses Matplotlib’s Computer Modern lookalike
    })

    amount = len(interested_fields) + (1 if include_cfl else 0) + (1 if include_dt else 0)
    fig, ax = plt.subplots(amount, 1, figsize=(figsize[0], figsize[1] * amount), sharex=True)
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
            ax[current_ax].plot(frames, linewidth=2.5, label=f'{key}{data_name}', zorder=2)
            if include_extra_stats:
                ax[current_ax].axhline(mean, color=COLOR_THEMES[i % len(COLOR_THEMES)]["mean"], linestyle='--', linewidth=1.2, label=f'Mean{data_name}: {mean:.2f}', zorder=1)
                ax[current_ax].axhline(median, color=COLOR_THEMES[i % len(COLOR_THEMES)]["median"], linestyle='--', linewidth=1.2, label=f'Median{data_name}: {median:.2f}', zorder=1)
                ax[current_ax].axhline(minimum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["min"], linestyle=':', linewidth=1.2, label=f'Min{data_name}: {minimum:.2f}', zorder=1)
                ax[current_ax].axhline(maximum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["max"], linestyle=':', linewidth=1.2, label=f'Max{data_name}: {maximum:.2f}', zorder=1)
            ax[current_ax].set_ylabel(key)
            ax[current_ax].set_title(f"{key} Over Time")
            ax[current_ax].legend(loc='best')
            ax[current_ax].grid(True)
            ax[current_ax].ticklabel_format(style='plain', axis='y', useOffset=False)
            
        current_ax += 1

    plt.xlabel("Timesteps")
    if include_title:
        fig.suptitle(title, fontsize="x-large")
        plt.tight_layout(rect=[0, 0, 1, .99])
    else:
        plt.tight_layout()
    plt.savefig(export_path, bbox_inches="tight")

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


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def create_combined_graph2(data_array, data_names, interested_fields, title,
                          include_title=True, include_cfl=True, include_dt=True,
                          include_extra_stats=True, export_path="./exports/combined.png",
                          figsize=(12,4), outlier_radius=100):
    """
    Creates and saves a combined multi-line graph showing the evolution of various fields over time.
    Automatically uses a broken y-axis if a dataset deviates too much from its initial value.

    Parameters:
        data_array (list of dict): List of datasets, each containing "frame_data" with field values per timestep.
        data_names (list of str): Names corresponding to each dataset for labeling in the graph legend.
        interested_fields (list of tuples): List of (field_name, target_key (like sum, max, ...)) pairs to extract from each frame.
        title (str): The title for the entire figure.
        include_cfl (bool): If True, includes the CFL field in the plots.
        include_dt (bool): If True, includes the timestep (dt) in the plots.
        include_extra_stats (bool): If True, overlays mean, median, min, and max lines for each dataset.
        export_path (str): File path where the final plot image will be saved.
        figsize (tuple): Size of one individual graph.
        outlier_radius (float): Threshold above or below first datapoint to consider a dataset an outlier.
    """
    frames_sets = []

    for data in data_array:
        frames = data["frame_data"]
        if not frames:
            print("No frame data collected")
            return
        
        field_frames = {}
        if include_cfl:
            field_frames["cfl"] = []
        if include_dt:
            field_frames["dt"] = []
        for interested_field, target in interested_fields:
            field_frames[f"{interested_field} {target}"] = []
              
        for key, grid in frames.items():
            if include_cfl:
                field_frames["cfl"].append(grid["cfl"])
            
            if include_dt:
                field_frames["dt"].append(grid["dt"])

            for interested_field, target in interested_fields:
                field_frames[f"{interested_field} {target}"].append(grid[interested_field][target])
        
        frames_sets.append(field_frames)

    plt.style.use("seaborn-v0_8-notebook")
    mpl.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "mathtext.fontset": "cm"
    })

    amount = len(interested_fields) + (1 if include_cfl else 0) + (1 if include_dt else 0)
    fig, ax = plt.subplots(amount, 1, figsize=(figsize[0], figsize[1] * amount), sharex=True)
    ax = np.atleast_1d(ax)

    current_ax = 0
    for key, _ in frames_sets[0].items():
        # Detect outliers based on deviation from first value
        outlier_idx = []
        for i in range(len(data_array)):
            frames = np.array(frames_sets[i][key])
            start = frames[0]
            if np.max(frames) > start + outlier_radius or np.min(frames) < start - outlier_radius:
                outlier_idx.append(i)

        # Determine axes
        if outlier_idx:
            # Broken y-axis
            fig2, (ax_top, ax_bottom) = plt.subplots(2,1,sharex=True,
                                                     figsize=(figsize[0], figsize[1]*amount),
                                                     gridspec_kw={'height_ratios':[1,3]})
            ax_bottom = np.atleast_1d(ax_bottom)
            ax_top = np.atleast_1d(ax_top)
            
            # Determine normal vs outlier range
            normal_vals = [np.max(frames_sets[i][key]) for i in range(len(data_array)) if i not in outlier_idx]
            normal_max = max(normal_vals) if normal_vals else 0
            normal_min = min(normal_vals) if normal_vals else 0
            top_vals = [np.max(frames_sets[i][key]) for i in outlier_idx]
            bottom_vals = [np.min(frames_sets[i][key]) for i in outlier_idx]
            top_min = min(top_vals) if top_vals else normal_max
            top_max = max(top_vals) if top_vals else normal_max
            bottom_min = min(normal_min, min(bottom_vals) if bottom_vals else normal_min)
            bottom_max = normal_max

            ax_bottom[0].set_ylim(bottom_min*0.95, bottom_max*1.05)
            ax_top[0].set_ylim(top_min*0.95, top_max*1.05)

            for i in range(len(data_array)):
                frames = np.array(frames_sets[i][key])
                label_name = f", {data_names[i]}" if data_names[i] != "" else ""
                if i in outlier_idx:
                    ax_top[0].plot(frames, linestyle=":", linewidth=2.5, label=f"{key}{label_name}", zorder=2)
                    ax_bottom[0].plot(frames, linestyle=":", linewidth=2.5, label=f"{key}{label_name}", zorder=2)
                else:
                    ax_bottom[0].plot(frames, linewidth=2.5, label=f"{key}{label_name}", zorder=2)

            ax_top[0].spines['bottom'].set_visible(False)
            ax_bottom[0].spines['top'].set_visible(False)
            ax_top[0].tick_params(labelbottom=False)
            ax_bottom[0].xaxis.tick_bottom()
            ax_bottom[0].set_ylabel(key)
            ax_top[0].set_title(f"{key} Over Time")
            ax_top[0].legend(loc='best')
            ax_bottom[0].grid(True)
            ax_bottom[0].ticklabel_format(style='plain', axis='y', useOffset=False)

        else:
            # Normal single axis
            for i in range(len(data_array)):
                frames = np.array(frames_sets[i][key])
                label_name = f", {data_names[i]}" if data_names[i] != "" else ""
                ax[current_ax].plot(frames, linewidth=2.5, label=f"{key}{label_name}", zorder=2)
            ax[current_ax].set_ylabel(key)
            ax[current_ax].set_title(f"{key} Over Time")
            ax[current_ax].legend(loc='best')
            ax[current_ax].grid(True)
            ax[current_ax].ticklabel_format(style='plain', axis='y', useOffset=False)

        current_ax += 1

    plt.xlabel("Timesteps")
    if include_title:
        fig.suptitle(title, fontsize="x-large")
        plt.tight_layout(rect=[0, 0, 1, .99])
    else:
        plt.tight_layout()
    plt.savefig(export_path, bbox_inches="tight")
