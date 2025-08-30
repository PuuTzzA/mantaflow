from pathlib import Path
import sys
import json
from graph_creation import *


BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()
FILENAME = "data.json"
DASHED = (0, (11, 11))

data_zalesak_rotation = [
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial_no_clamped_redistro',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear',
    '/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial',

    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic_no_clamped_redistro',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic',
    #'/2_fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite',
]

data_shear_flow = [
    '/1_fixed_vel_shear_flow/shear_flow_traditional_RK4_2_polynomial',
    '/1_fixed_vel_shear_flow/shear_flow_conserving_2_polynomial_no_clamped_redistro',
    '/1_fixed_vel_shear_flow/shear_flow_conserving_0_linear',
    '/1_fixed_vel_shear_flow/shear_flow_conserving_2_polynomial',
    '/1_fixed_vel_shear_flow/shear_flow_traditional_RK4_0_linear',

    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic_no_clamped_redistro',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_conserving_1_cubic',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE1',
    #'/1_fixed_vel_shear_flow_low_cfl/shear_flow_traditional_EE2',
    #'/1_fixed_vel_shear_flow/shear_flow_traditional_RK4_1_cubic',
    #'/1_fixed_vel_shear_flow/shear_flow_traditional_RK4_3_monotone_hermite'
]

data_simple_plume_cfl_5 = [
    '/3_simple_plume/cfl_05/plume_2d_low_traditional_RK4_2_polynomial',
    '/3_simple_plume/cfl_05/plume_2d_low_conserving_2_polynomial_no_clamped_redistro',
    '/3_simple_plume/cfl_05/plume_2d_low_traditional_RK4_0_linear',

    '/3_simple_plume/cfl_05/plume_2d_low_conserving_0_linear',
    '/3_simple_plume/cfl_05/plume_2d_low_conserving_2_polynomial',
]

data_simple_obstacle_cfl_5 = [
    '/4_simple_obstacle/cfl_05/obstacle_2d_low_traditional_RK4_2_polynomial',
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

data_different_cfl = [
    '/6_different_cfl/40_dif_cfl_conserving_polynomial_local_cfl',
    '/6_different_cfl/40_dif_cfl_conserving_polynomial',
    '/6_different_cfl/40_dif_cfl_traditional_polynomial_local_cfl',
    '/6_different_cfl/40_dif_cfl_traditional_polynomial',

    '/6_different_cfl/30_dif_cfl_conserving_polynomial_local_cfl',
    '/6_different_cfl/30_dif_cfl_conserving_polynomial',
    '/6_different_cfl/30_dif_cfl_traditional_polynomial_local_cfl',
    '/6_different_cfl/30_dif_cfl_traditional_polynomial',

    '/6_different_cfl/20_dif_cfl_conserving_polynomial_local_cfl',
    '/6_different_cfl/20_dif_cfl_conserving_polynomial',
    '/6_different_cfl/20_dif_cfl_traditional_polynomial_local_cfl',
    '/6_different_cfl/20_dif_cfl_traditional_polynomial',

    '/6_different_cfl/10_dif_cfl_conserving_polynomial_local_cfl',
    '/6_different_cfl/10_dif_cfl_conserving_polynomial',
    '/6_different_cfl/10_dif_cfl_traditional_polynomial_local_cfl',
    '/6_different_cfl/10_dif_cfl_traditional_polynomial',

    '/6_different_cfl/05_dif_cfl_conserving_polynomial_local_cfl',
    '/6_different_cfl/05_dif_cfl_conserving_polynomial',
    '/6_different_cfl/05_dif_cfl_traditional_polynomial_local_cfl',
    '/6_different_cfl/05_dif_cfl_traditional_polynomial',

#    '/6_different_cfl/01_dif_cfl_conserving_polynomial_local_cfl',
#    '/6_different_cfl/01_dif_cfl_conserving_polynomial',
#    '/6_different_cfl/01_dif_cfl_traditional_polynomial_local_cfl',
#    '/6_different_cfl/01_dif_cfl_traditional_polynomial',
]

data = data_shear_flow
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


data_2 = data_zalesak_rotation
input_datas_2 = []
for d in data_2:
    path = (BASEDIR / d.strip('/') / FILENAME).resolve()
    with open(path) as f:
	    input_datas_2.append(json.load(f))



""" #1_shear_flow
input_titles = ["not conserving, polynomial interpolation", "conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "not conserving, linear interpolation"]
title = "shear_flow_total_mass_2"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [4, 0, 2, 1, 3]
#labelOrder = None

linestyles = ['solid', DASHED, 'solid', DASHED, 'solid']
linewidths = [5.5]
blue =      plt.cm.tab10.colors[0]
orange =    plt.cm.tab10.colors[1] 
green =     plt.cm.tab10.colors[2]
red =       plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [green, red, blue, orange, violett]
#colors = plt.cm.tab10.colors

margins = (0.1, 0.1)
figsize = (11.8, 8.5)
allTextSize = 24

bBoxAnchor = (0.5, 1.3)
extraMargins = (0, 0, 0.05, 0.3) #left, right, bottom, top

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins,
                      allTextSize=allTextSize, bBoxAnchor=bBoxAnchor, extraMargins=extraMargins, data_2_share_yAxis=False)
 """

""" #2_zalesak_rotation
interested_fields = [["testField", "sum"]]
#input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "shear_flow_total_mass2"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [3, 1, 0, 2]
labelOrder = None

linestyles = ['solid', 'solid', DASHED, 'solid', DASHED]
linewidths = [5]
blue =  plt.cm.tab10.colors[0]
orange =plt.cm.tab10.colors[1] 
green = plt.cm.tab10.colors[2]
red =   plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [violett, green, red, blue, orange]

margins = (0.086666, 0.13)
figsize = (11.8, 7.5)
allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, allTextSize=allTextSize)
 """

#shear flow and zalesak rotation
input_titles = ["not conserving, polynomial interpolation", "conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "not conserving, linear interpolation"]
title = "shear_flow_zalesak_total_mass"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [4, 0, 2, 1, 3]
#labelOrder = None

linestyles = ['solid', DASHED, 'solid', DASHED, 'solid']
linewidths = [5.5]
blue =      plt.cm.tab10.colors[0]
orange =    plt.cm.tab10.colors[1] 
green =     plt.cm.tab10.colors[2]
red =       plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [green, red, blue, orange, violett]
#colors = plt.cm.tab10.colors

linestyles2 = ['solid', 'solid', DASHED, 'solid', DASHED]
colors2 = [violett, green, red, blue, orange]
linewidths2 = [5.5]

figsize = (7.7, 6)
allTextSize = 24

bBoxAnchor = (0.525, 1.207)
margins = (0.1, 0.1, 0.128, 0.6) #left, right, bottom, top
margins2 = (0.1, 0.1, 0.128, 0.4)

create_combined_graph(data_array=input_datas, data_array_2=input_datas_2, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Relative change in mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors,
                      allTextSize=allTextSize, bBoxAnchor=bBoxAnchor, margins=margins, data_2_share_yAxis=False, 
                      linestyles2=linestyles2, colors2=colors2, linewidths2=linewidths2, margins2=margins2)


""" #3_simple_plume_cfl_5
#input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "simple_plume_cfl_05"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [1, 2, 0, 3]
labelOrder = None

linestyles = ['solid', DASHED, 'solid', 'solid', DASHED]
linewidths = [5.5]
blue =      plt.cm.tab10.colors[0]
orange =    plt.cm.tab10.colors[1] 
green =     plt.cm.tab10.colors[2]
red =       plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [green, red, violett, blue, orange]

margins = (0.1, 0.1, 0.15, 0.15) #left, right, bottom, top
figsize = (7.5, 5)

allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, 
                      allTextSize=allTextSize, data_array_2=None)
 """

""" #3_simple_obstacle_cfl_5
#input_titles = ["conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "semi-Lagrangian, linear interpolation"]
title = "simple_obstacle_cfl_05"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [1, 2, 0, 3]
labelOrder = None

linestyles = ['solid', DASHED, 'solid', 'solid', DASHED]
linewidths = [5.5]
blue =      plt.cm.tab10.colors[0]
orange =    plt.cm.tab10.colors[1] 
green =     plt.cm.tab10.colors[2]
red =       plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [green, red, violett, blue, orange]

margins = (0.1, 0.1, 0.15, 0.15) #left, right, bottom, top
figsize = (7.5, 5)

allTextSize = 24

create_combined_graph(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors, margins=margins, 
                      allTextSize=allTextSize, data_array_2=None)
 """

""" #simple plume and simple obstacle
input_titles = ["not conserving, polynomial interpolation", "conserving, polynomial interpolation, no redistribution", "not conserving, linear interpolation", "conserving, linear interpolation", "conserving, polynomial interpolation"]
title = "simple_plume_and_obstacle_total_mass"
output_path = input_paths[0].parent.parent / f"{title}.pdf"
labelOrder = [2, 0, 3, 1, 4]
#labelOrder = None

linestyles = ['solid', DASHED, 'solid', 'solid', DASHED]
linewidths = [5.5]
blue =      plt.cm.tab10.colors[0]
orange =    plt.cm.tab10.colors[1] 
green =     plt.cm.tab10.colors[2]
red =       plt.cm.tab10.colors[3]
violett =   plt.cm.tab10.colors[4] 
colors = [green, red, violett, blue, orange]
#colors = plt.cm.tab10.colors

figsize = (7.7, 6)
allTextSize = 24

bBoxAnchor = (0.525, 1.207)
margins = (0.1, 0.1, 0.128, 0.6) #left, right, bottom, top
margins2 = (0.1, 0.1, 0.128, 0.6)

create_combined_graph(data_array=input_datas, data_array_2=input_datas_2, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, figsize=figsize, 
                      yAxisLabel="Total mass", labelOrder=labelOrder, linestyles=linestyles, linewidths=linewidths, colors=colors,
                      allTextSize=allTextSize, bBoxAnchor=bBoxAnchor, margins=margins, data_2_share_yAxis=True, 
                      linestyles2=linestyles, colors2=colors, linewidths2=linewidths, margins2=margins2)
 """


""" # different cfl 2d
#input_titles = ["not conserving, polynomial interpolation", "conserving, polynomial interpolation, no redistribution", "conserving, linear interpolation", "conserving, polynomial interpolation", "not conserving, linear interpolation"]
title = "different_cfl"
output_path = input_paths[0].parent.parent / f"{title}.pdf"

create_combined_graph_old(data_array=input_datas, data_names=input_titles, interested_fields=interested_fields, 
                      title=title, include_title=False, include_cfl=False, include_dt=False, include_extra_stats=False, export_path=output_path, 
                      figsize=(12, 6))
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