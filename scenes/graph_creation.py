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


def create_combined_graph_old(data_array, data_names, interested_fields, title, include_title=True, include_cfl=True, include_dt=True, include_extra_stats=True, export_path="./exports/combinded.png", figsize=(12,4)):
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


def create_combined_graph(data_array, data_names, interested_fields, title, include_title=True, include_cfl=True, include_dt=True, include_extra_stats=True, export_path="./exports/combinded.png", figsize=(12,4), yAxisLabel = "---", 
                          labelOrder = None, linestyles=['solid'], linewidths=[3.5], colors=plt.cm.tab10.colors[0], margins=(.1, .1)):
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
        field_frames["time"] = []
        if include_cfl:
            field_frames["cfl"] = []
        if include_dt:
            field_frames["dt"] = []
        for interested_field, target in interested_fields:
            field_frames[f"{interested_field} {target}"] = []
              
        iter = 0
        for key, grid in frames.items():
            field_frames["time"].append(float(grid["dt"]) + (0 if iter == 0 else field_frames["time"][iter - 1]))
            iter += 1

            if include_cfl:
                field_frames["cfl"].append(grid["cfl"])
            
            if include_dt:
                field_frames["dt"].append(grid["dt"])

            for interested_field, target in interested_fields:
                field_frames[f"{interested_field} {target}"].append(grid[interested_field][target])
        
        frames_sets.append(field_frames)

    #plt.style.use("seaborn-v0_8-deep")
    #plt.style.use("seaborn-v0_8-notebook")
    plt.style.use("classic")

    plt.rcParams["legend.fontsize"] = "medium"  # or a number like 10
    plt.rcParams["axes.labelsize"] = "medium"

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
        if key == "time":
            continue

        for i in range(len(data_array)):     

            frames = np.array(frames_sets[i][key])
            minimum = np.min(frames)
            maximum = np.max(frames)
            mean = np.mean(frames)
            std_dev = np.std(frames)
            median = np.median(frames)

            data_name = data_names[i]
            ax[current_ax].plot(np.array(frames_sets[i]["time"]), frames, 
                                color=colors[i % len(colors)], linewidth=linewidths[i % len(linewidths)], linestyle=linestyles[i % len(linestyles)], label=f'{data_name}', zorder=2)
            if include_extra_stats:
                ax[current_ax].axhline(mean, color=COLOR_THEMES[i % len(COLOR_THEMES)]["mean"], linestyle='--', linewidth=1.2, label=f'Mean{data_name}: {mean:.2f}', zorder=1)
                ax[current_ax].axhline(median, color=COLOR_THEMES[i % len(COLOR_THEMES)]["median"], linestyle='--', linewidth=1.2, label=f'Median{data_name}: {median:.2f}', zorder=1)
                ax[current_ax].axhline(minimum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["min"], linestyle=':', linewidth=1.2, label=f'Min{data_name}: {minimum:.2f}', zorder=1)
                ax[current_ax].axhline(maximum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["max"], linestyle=':', linewidth=1.2, label=f'Max{data_name}: {maximum:.2f}', zorder=1)
            ax[current_ax].set_ylabel(yAxisLabel)
            #ax[current_ax].set_title(f"{key} Over Time")

            ax[current_ax].set_axisbelow(True)
            ax[current_ax].grid(True, linestyle="-", color='gray', linewidth=0.5, alpha=0.3)
            ax[current_ax].ticklabel_format(style='plain', axis='y', useOffset=False)
            
        current_ax += 1
      
    handles, labels = plt.gca().get_legend_handles_labels()

    if labelOrder is not None:
        plt.legend([handles[i] for i in labelOrder], [labels[i] for i in labelOrder], framealpha=0.5, edgecolor='black', loc="upper left", bbox_to_anchor=(1.05, 1))
    else:
        plt.legend(loc='best')

    plt.margins(x=margins[0], y=margins[1])  # 5% padding on both axes

    plt.xlabel("Time")
    if include_title:
        fig.suptitle(title, fontsize="x-large")
        plt.tight_layout(rect=[0, 0, 1, .99])
    else:
        plt.tight_layout()
    plt.savefig(export_path, bbox_inches="tight")


def create_combined_graph_with_clipping(
    data_array, data_names, interested_fields, title,
    include_title=True, include_cfl=True, include_dt=True, include_extra_stats=True,
    export_path="./exports/combinded.png", figsize=(12,4),
    r1=None, r2=None
):
    """
    Creates and saves a combined multi-line graph showing the evolution of various fields over time.
    Supports clipping and dashed line rendering for outliers.

    Parameters:
        ...
        r1 (float): Threshold for dashed line (relative to first value).
        r2 (float): Threshold for clipping (relative to first value).
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
        for i in range(len(data_array)):     
            frames = np.array(frames_sets[i][key])
            x_vals = np.arange(len(frames))
            first_val = frames[0]

            # Segment the data into normal / dashed / clipped
            normal_x, normal_y = [], []
            dashed_x, dashed_y = [], []
            clipped = False

            for t, val in zip(x_vals, frames):
                diff = abs(val - first_val)

                if r2 is not None and diff > r2:
                    # Stop plotting completely
                    clipped = True
                    break
                elif r1 is not None and diff > r1:
                    dashed_x.append(t)
                    dashed_y.append(val)
                else:
                    normal_x.append(t)
                    normal_y.append(val)

            data_name = f", {data_names[i]}" if data_names[i] != "" else ""

            # Plot normal part
            if normal_x:
                ax[current_ax].plot(normal_x, normal_y, linewidth=2.5,
                                    label=f'{key}{data_name}', zorder=2)

            # Plot dashed part
            if dashed_x:
                ax[current_ax].plot(dashed_x, dashed_y, linewidth=2.5,
                                    linestyle="--", label=f'{key}{data_name} (outlier)',
                                    zorder=2)

            # Extra stats (mean/median/min/max) only if not clipped
            if include_extra_stats and not clipped:
                mean = np.mean(frames)
                median = np.median(frames)
                minimum = np.min(frames)
                maximum = np.max(frames)

                ax[current_ax].axhline(mean, color=COLOR_THEMES[i % len(COLOR_THEMES)]["mean"],
                                       linestyle='--', linewidth=1.2, label=f'Mean{data_name}: {mean:.2f}', zorder=1)
                ax[current_ax].axhline(median, color=COLOR_THEMES[i % len(COLOR_THEMES)]["median"],
                                       linestyle='--', linewidth=1.2, label=f'Median{data_name}: {median:.2f}', zorder=1)
                ax[current_ax].axhline(minimum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["min"],
                                       linestyle=':', linewidth=1.2, label=f'Min{data_name}: {minimum:.2f}', zorder=1)
                ax[current_ax].axhline(maximum, color=COLOR_THEMES[i % len(COLOR_THEMES)]["max"],
                                       linestyle=':', linewidth=1.2, label=f'Max{data_name}: {maximum:.2f}', zorder=1)

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

def create_combined_graph_broken_axis(
    data_array, data_names, interested_fields, title,
    include_title=True, include_cfl=True, include_dt=True, include_extra_stats=True,
    export_path="./exports/combinded.png", figsize=(12,4),
    y_ranges=None, y_gap=0.05, smallTitles = ["moin"], yLabels = ["meister"]
):
    """
    Creates and saves a combined multi-line graph with broken y-axis support.

    Parameters:
        ...
        y_ranges (list of [low, high]): Visible y-intervals for broken axis.
        figsize (tuple or list of tuples): If tuple, applies same size to all.
                                           If list, individual sizes per y_range (height ratios).
        y_gap (float): Vertical gap between broken y-axis parts (like hspace). 0 means no gap.
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

    if isinstance(figsize, list):
        width = figsize[0][0]
        heights = [h for (_,h) in figsize]
        height_ratios = heights
        fig_height = sum(heights) * amount
        fig, axes = plt.subplots(amount * len(y_ranges), 1,
                                 figsize=(width, fig_height),
                                 sharex=True,
                                 gridspec_kw={'height_ratios': height_ratios * amount,
                                              'hspace': y_gap})
    else:
        fig, axes = plt.subplots(amount * (len(y_ranges) if y_ranges else 1), 1,
                                 figsize=(figsize[0], figsize[1] * amount * (len(y_ranges) if y_ranges else 1)),
                                 sharex=True,
                                 gridspec_kw={'hspace': y_gap})
    axes = np.atleast_1d(axes)

    current_ax = 0
    for key, _ in frames_sets[0].items():
        handles, labels = [], []

        for i in range(len(data_array)):     
            frames = np.array(frames_sets[i][key])
            x_vals = np.arange(len(frames))
            data_name = f"{data_names[i]}" if data_names[i] != "" else ""

            if y_ranges is None:
                line, = axes[current_ax].plot(x_vals, frames, linewidth=2.5,
                                              label=f'{data_name}', zorder=2)
                handles.append(line)
                labels.append(f'{data_name}')
            else:
                for j, (low, high) in enumerate(y_ranges):
                    ax = axes[current_ax * len(y_ranges) + j]
                    line, = ax.plot(x_vals, frames, linewidth=2.5,
                                    label=f'{data_name}', zorder=2)
                    ax.set_ylim(low, high)
                    handles.append(line)
                    labels.append(f'{data_name}')

        # Deduplicate legend entries
        unique = {}
        for h, l in zip(handles, labels):
            unique[l] = h
        handles, labels = list(unique.values()), list(unique.keys())

        if y_ranges is None:
            ax = axes[current_ax]
            ax.set_ylabel(f"{yLabels[0]}")
            ax.set_title(f"{smallTitles[0]}")
            ax.grid(True)
            ax.ticklabel_format(style='plain', axis='y', useOffset=False)
        else:
            # Calculate the vertical center position for this group of subplots
            group_start = current_ax * len(y_ranges)
            group_end = (current_ax + 1) * len(y_ranges) - 1
            
            # Get the positions of the first and last axes in this group
            ax_top = axes[group_start]
            ax_bottom = axes[group_end]
            
            # Calculate center position in figure coordinates
            top_pos = ax_top.get_position().y1
            bottom_pos = ax_bottom.get_position().y0
            center_y = (top_pos + bottom_pos) / 2
            
            # Add centered y-label for this group
            fig.text(0.04, center_y, f"{yLabels[current_ax] if current_ax < len(yLabels) else yLabels[0]}", 
                    rotation=90, va='center', ha='center', fontsize=plt.rcParams['axes.labelsize'])
            
            for j, (low, high) in enumerate(y_ranges):
                ax = axes[current_ax * len(y_ranges) + j]
                # Don't set individual y-labels anymore
                ax.grid(True)
                ax.ticklabel_format(style='plain', axis='y', useOffset=False)

            # Hide spines & add diagonal cut marks
            for j in range(len(y_ranges)-1):
                ax_top = axes[current_ax * len(y_ranges) + j]
                ax_bottom = axes[current_ax * len(y_ranges) + j + 1]
                ax_top.spines.bottom.set_visible(False)
                ax_bottom.spines.top.set_visible(False)
                ax_top.xaxis.tick_top()
                ax_top.tick_params(labeltop=False)
                ax_bottom.xaxis.tick_bottom()

                d = .5
                kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
                              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
                ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
                ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)

        current_ax += 1

    plt.xlabel("Step")
    
    # Collect all handles & labels from every axis
    all_handles, all_labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        all_handles.extend(h)
        all_labels.extend(l)

    # Deduplicate
    unique = {}
    for h, l in zip(all_handles, all_labels):
        unique[l] = h
    handles, labels = list(unique.values()), list(unique.keys())

    # Place legend inside the plot area on the first axis
    if handles and labels:
        # Use the first (top) axis for legend placement
        legend_ax = axes[0]
        
        # Simple legend without any styling
        legend_ax.legend(handles, labels, loc='best')

    if include_title:
        fig.suptitle(title, fontsize="x-large")
        if y_ranges is None:
            plt.tight_layout(rect=[0, 0, 1, 0.96])
        else:
            plt.subplots_adjust(top=0.92, left=0.1)  # Add left margin for y-labels
    else:
        if y_ranges is None:
            plt.tight_layout()
        else:
            plt.subplots_adjust(left=0.1)  # Add left margin for y-labels

    # Save the figure
    plt.savefig(export_path, bbox_inches="tight", dpi=300)
    plt.close(fig)  # Close the figure to free memory