import subprocess
from pathlib import Path

# Constants
BASEDIR = (Path(__file__).parent.parent / "exports/").resolve()

videos_zalesak_rotation = [
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_0_linear/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_conserving_2_polynomial/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE1/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_EE2/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_0_linear/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_1_cubic/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_2_polynomial/',
    '/fixed_vel_zalesak_rotation/zalesak_rotation_traditional_RK4_3_monotone_hermite/',
]

videos_shear_flow = [
    '/fixed_vel_shear_flow/shear_flow_conserving_0_linear',
    '/fixed_vel_shear_flow/shear_flow_conserving_1_cubic',
    '/fixed_vel_shear_flow/shear_flow_conserving_2_polynomial',
    '/fixed_vel_shear_flow/shear_flow_traditional_EE1',
    '/fixed_vel_shear_flow/shear_flow_traditional_EE2',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_0_linear',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_1_cubic',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_2_polynomial',
    '/fixed_vel_shear_flow/shear_flow_traditional_RK4_3_monotone_hermite',
]

videos_plume_2d_low = [
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_0_linear',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_1_cubic',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_conserving_2_polynomial',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE1',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_EE2',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_0_linear',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_1_cubic',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_2_polynomial',
    '/simple_plume_obstacle_2d_low/simple_plume_obstacle_2d_low_traditional_RK4_3_monotone_hermite',
]

videos_plume_2d_high = [
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_0_linear',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_1_cubic',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_conserving_2_polynomial',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE1',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_EE2',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_0_linear',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_1_cubic',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_2_polynomial',
    '/simple_plume_obstacle_2d_high/simple_plume_obstacle_2d_high_traditional_RK4_3_monotone_hermite',
]

videos = videos_plume_2d_high
FILENAME = "testField.webm" # for fixed vel
FILENAME = "density.webm" # for plume

# Resolve full paths
input_paths = []
for v in videos:
    path = (BASEDIR / v.strip('/') / FILENAME).resolve()
    input_paths.append(path)

# Construct xstack layout
layout = (
    "0_0|w0_0|w0+w1_0|"
    "0_h0|w3_h0|w3+w4_h0|"
    "0_h0+h3|w6_h0+h3|w6+w7_h0+h3"
)

# Build ffmpeg xstack command
xstack_inputs = []
for p in input_paths:
    xstack_inputs += ["-i", str(p)]

filter_complex = f"xstack=inputs=9:layout={layout}"

output_path = input_paths[0].parent.parent / f"combined_{FILENAME.removesuffix(".webm")}.webm"

cmd = [
    "ffmpeg",
    "-y",  # Overwrite
    *xstack_inputs,
    "-filter_complex", filter_complex,
    "-c:v", "libvpx-vp9",
    "-crf", "30",
    "-b:v", "0",
    str(output_path)
]

print("Creating 3x3 grid...")
subprocess.run(cmd, check=True)

print(f"Grid video saved to: {output_path}")
