from pathlib import Path
import sys
import json
from graph_creation import *


BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
FILENAME = "data.json"

data_zalesak_rotation = [
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial_no_clamped_redistro',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear',

    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic_no_clamped_redistro',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic.json',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite',
]

data_shear_flow_low = [
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_2_polynomial_no_clamped_redistro',
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_0_linear',
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_2_polynomial',
    '/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_0_linear',

    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic_no_clamped_redistro',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE1',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE2',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_1_cubic',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_2_polynomial',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_RK4_3_monotone_hermitn'
]

data_simple_plume_cfl_5 = [
    '/3_simple_plume/cfl_05/plume_2d_low_conserving_2_polynomial_no_clamped_redistro',
    '/3_simple_plume/cfl_05/plume_2d_low_traditional_RK4_0_linear',

    '/3_simple_plume/cfl_05/plume_2d_low_conserving_0_linear',
    '/3_simple_plume/cfl_05/plume_2d_low_conserving_2_polynomial',
]

data_simple_obstacle_cfl_5 = [
    '/4_simple_obstacle/cfl_05/obstacle_2d_low_conserving_2_polynomial_no_clamped_redistro',
    '/4_simple_obstacle/cfl_05/obstacle_2d_low_traditional_RK4_0_linear',
    '/4_simple_obstacle/cfl_05/obstacle_2d_low_conserving_0_linear',
    '/4_simple_obstacle/cfl_05/obstacle_2d_low_conserving_2_polynomial',
]

data_simple_plume_cfl_20 = [
    '/3_simple_plume/cfl_20/plume_2d_high_conserving_0_linear',
    '/3_simple_plume/cfl_20/plume_2d_high_conserving_2_polynomial',
    '/3_simple_plume/cfl_20/plume_2d_high_conserving_2_polynomial_local_cfl',
    '/3_simple_plume/cfl_20/plume_2d_high_traditional_RK4_0_linear',
]

data_simple_obstacle_cfl_30 = [
    '/4_simple_obstacle/cfl_30/obstacle_2d_high_conserving_0_linear',
    '/4_simple_obstacle/cfl_30/obstacle_2d_high_conserving_2_polynomial',
    '/4_simple_obstacle/cfl_30/obstacle_2d_high_conserving_2_polynomial_local_cfl',
    '/4_simple_obstacle/cfl_30/obstacle_2d_high_traditional_RK4_0_linear',
]

data = data_shear_flow_low
interested_fields = [["testField", "sum"]] # for fixed vel Field
#interested_fields = [["fixed_volume", "sum"]] # for plume

# Resolve full paths
input_paths = []
input_datas = []
input_titles = []
for d in data:
    path = (BASEDIR / d.strip('/') / FILENAME).resolve()
    title = path.parent.name

    input_titles.append(title)
    input_paths.append(path)

    with open(path) as f:
	    input_datas.append(json.load(f))


data_2 = data_shear_flow_low
input_datas_2 = []
for d in data_2:
    path = (BASEDIR / d.strip('/') / FILENAME).resolve()
    with open(path) as f:
	    input_datas_2.append(json.load(f))


# Create combined Graph

#1_shear_flow
input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "shear_flow_total_mass_2"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [3, 1, 0, 2]

linestyles = ['solid', 'solid', (0, (7.75, 7.75)), 'solid']
linewidths = [5.5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
colors = [green, blue, orange, red]

margins = (0.1, 0.1)
figsize = (11.8, 7.5)
allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, allTextSize=allTextSize)



""" #2_zalesak_rotation
interested_fields = [["testField", "sum"]]
input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "shear_flow_total_mass2"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [3, 1, 0, 2]

linestyles = ['solid', 'solid', 'dashed', 'solid']
linewidths = [5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
colors = [green, blue, orange, red]

margins = (0.086666, 0.13)
figsize = (7.5, 5)

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins)
 """


""" #3_simple_plume_cfl_5
#input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "simple_plume_cfl_05"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [1, 2, 0, 3]
#labelOrder = None

linestyles = ['solid', 'solid', 'solid', (0, (7.75, 7.75))]
linewidths = [5.5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
colors = [green, red, blue, orange]

margins = (0.1, 0.15)
figsize = (7.5, 5)

allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, 
                      allTextSize=allTextSize, data_array_2=None)
 """


""" #3_simple_obstacle_cfl_5
input_titles = ["conserving, polynomial interpolation, no redistribution", "semi-Lagrangian, linear interpolation", "conserving, linear interpolation", "conserving, polynomial interpolation"]
title = "legend_2"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [1, 2, 0, 3]
#labelOrder = None

linestyles = ['solid', 'solid', 'solid', (0, (7.75, 7.75))]
linewidths = [5.5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
colors = [green, red, blue, orange]

margins = (0.1, 0.15)
figsize = (7.5, 5)

allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, allTextSize=allTextSize)
 """


""" # plume and obstacle high
#input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "simple_obstacle_cfl_30"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [1, 2, 0, 3]
#labelOrder = None

linestyles = ['solid', 'solid', 'solid', (0, (7.75, 5.5))]
linewidths = [5.5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
colors = [green, red, blue, orange]

margins = (0.086666, 0.13)
figsize = (7.5, 5)

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins)
 """






























#Setup for Shear Flow
""" labels = ["Conserving, linear", "Conserving, polynomial", "EE, linear", "RK4, polynomial"]
output_path = input_paths[0].parent.parent / f"shear_flow_mass_plot.pdf"

create_combined_graph_broken_axis(data_array=input_datas, data_names=labels, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, 
                      figsize=[(12, 3), (12, 1)], y_ranges=[[1960, 2200], [0, 800]], y_gap=0.05, smallTitles=["Total Mass"], yLabels=["Mass"])
 """


#Setup for Zalesak Rotation
""" labels = ["Conserving, linear", "Conserving, polynomial", "EE, linear", "RK4, polynomial"]
output_path = input_paths[0].parent.parent / f"zalesak_rotation_mass_plot.pdf"

create_combined_graph_broken_axis(data_array=input_datas, data_names=labels, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, 
                      figsize=[(12, 3), (12, 1)], y_ranges=[[4306, 4335], [2000, 2900]], y_gap=0.05, smallTitles=["Total Mass"], yLabels=["Mass"])
 """