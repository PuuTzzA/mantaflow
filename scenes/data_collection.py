import os
from manta import *
from pathlib import Path
import subprocess
import shutil
import json
import numpy as np
import matplotlib.pyplot as plt
import math

class Data_collectior:
    def __init__(self, title="no_title_specified", base_dir="../analysis/experiments/", params=None, export_data=True, export_images=False, export_videos=False, trackable_grid_names=[], tracked_grids_indeces=[], fixed_volume="fixed_volume"):
        self.title = title
        self.base_dir = Path(base_dir).expanduser().resolve()
        self.stats_dir = self.base_dir / f"{self.title}_stats"

        if params is not None:
            self.dim = params["dimension"]
            self.res = [params["resolutionX"], params["resolutionY"], params["resolutionZ"]]
            self.maxCfl = params["maxCFL"]
            self.max_time = params["max_time"]
            self.dt = params["dt"]
        
        self.export_data = export_data
        self.export_images = export_images
        self.export_videos = export_videos
        self.trackable_grids=trackable_grid_names
        self.tracked_grids_indeces = tracked_grids_indeces
        self.fixed_volume = fixed_volume

    def init(self):
        if self.export_data or self.export_images:
            shutil.rmtree(self.stats_dir, ignore_errors=True)
        self.stats_dir.mkdir(parents=True, exist_ok=True)
        
        self.data = {}
        self.data["title"] = self.title
        self.data["dim"] = self.dim
        self.data["res"] = self.res
        self.data["maxCfl"] = self.maxCfl
        self.data["max_time"] = self.max_time
        self.data["dt"] = self.dt
        self.data["frame_data"] = {}

        self.current_frame = 0
        
        if (self.export_images):
           for index in self.tracked_grids_indeces:
            (self.stats_dir / f"{self.trackable_grids[index][0]}_frames").mkdir(parents=True, exist_ok=True)

    def step(self, solver, flags, vel, gui=None, windowSize=[800, 800], camPos=[0, 0, -1.3]):
        #self.current_frame = math.floor(solver.timeTotal)
        self.data["frame_data"][str(self.current_frame).zfill(4)] = {}
        self.data["frame_data"][str(self.current_frame).zfill(4)]["cfl"] = vel.getMaxAbs() * solver.timestep
        self.data["frame_data"][str(self.current_frame).zfill(4)]["dt"] = solver.timestep

        if self.export_data:
            for i in range(len(self.trackable_grids)):
                if i in self.tracked_grids_indeces:    
                    name = self.trackable_grids[i][0]              
                    grid = self.trackable_grids[i][1]
                    self.data["frame_data"][str(self.current_frame).zfill(4)][name] = json.loads(realGridStats(grid=grid, flags=flags))

        if self.export_images and gui is not None:
            gui.windowSize(windowSize[0], windowSize[1])
            gui.setCamPos(camPos[0], camPos[1], camPos[2])

            for i in range(len(self.trackable_grids)):
                if i in self.tracked_grids_indeces:
                    name = self.trackable_grids[i][0]
                    gui.screenshot(str(self.stats_dir / f"{name}_frames" / f"{name}_{str(math.floor(self.current_frame)).zfill(4)}.png"))
                    
                gui.nextRealGrid()
                gui.update()
        
        self.current_frame += 1

    def computeStats(self):
        frames = self.data["frame_data"]
        if not frames:
            print("No frame data collected")
            return

        self.data["results"] = {}
        cfl_frames = []

        for f_key, grids in frames.items():
            try:
                current_cfl = grids["cfl"]
            except (KeyError):
                continue
            cfl_frames.append(current_cfl)
        
        cfl_frames = np.array(cfl_frames)

        cfl_minimum = np.min(cfl_frames)
        cfl_maximum = np.max(cfl_frames)
        cfl_mean = np.mean(cfl_frames)
        cfl_std_dev = np.std(cfl_frames)
        cfl_median = np.median(cfl_frames)

        self.data["results"]["cfl"] = {}
        self.data["results"]["cfl"]["minCfl"] = cfl_minimum
        self.data["results"]["cfl"]["maxCfl"] = cfl_maximum
        self.data["results"]["cfl"]["average"] = cfl_mean
        self.data["results"]["cfl"]["std_dev"] = cfl_std_dev
        self.data["results"]["cfl"]["median"] = cfl_median

        fixed_volume_frames = []
        for index in self.tracked_grids_indeces:
            grid_name = self.trackable_grids[index][0]
            min_sum   = float("inf")
            max_sum   = -float("inf")
            found_any = False

            for f_key, grids in frames.items():
                grid_stats = grids.get(grid_name)
                if not grid_stats:
                    continue   # grid wasn’t recorded this frame
                try:
                    s = float(grid_stats["sum"])
                except (KeyError, TypeError, ValueError):
                    continue   # malformed entry; skip

                found_any = True
                min_sum = min(min_sum, s)
                max_sum = max(max_sum, s)  

                if grid_name == self.fixed_volume:
                    fixed_volume_frames.append(s)

            if not found_any:
                continue

            self.data["results"][grid_name] = {}
            self.data["results"][grid_name]["minSum"] = min_sum
            self.data["results"][grid_name]["maxSum"] = max_sum

        if (self.export_data):
            # Convert to numpy arrays
            cfl_frames = np.array(cfl_frames)
            cfl_minimum = np.min(cfl_frames)
            cfl_maximum = np.max(cfl_frames)
            cfl_mean = np.mean(cfl_frames)
            cfl_std_dev = np.std(cfl_frames)
            cfl_median = np.median(cfl_frames)

            has_fixed_volume = len(fixed_volume_frames) > 0
            if has_fixed_volume:
                fixed_volume_frames = np.array(fixed_volume_frames)
                fixed_volume_minimum = np.min(fixed_volume_frames)
                fixed_volume_maximum = np.max(fixed_volume_frames)
                fixed_volume_mean = np.mean(fixed_volume_frames)
                fixed_volume_std_dev = np.std(fixed_volume_frames)
                fixed_volume_median = np.median(fixed_volume_frames)

            # Create subplots with shared x-axis
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

            # --- Plot CFL ---
            ax1.plot(cfl_frames, linestyle='-', linewidth=2.5, color='blue', label='CFL', zorder=2)
            ax1.axhline(cfl_mean, color='purple', linestyle='--', linewidth=1.2, label=f'Mean: {cfl_mean:.2f}', zorder=1)
            ax1.axhline(cfl_median, color='deeppink', linestyle='--', linewidth=1.2, label=f'Median: {cfl_median:.2f}', zorder=1)
            ax1.axhline(cfl_minimum, color='darkslategray', linestyle=':', linewidth=1.2, label=f'Min: {cfl_minimum:.2f}', zorder=1)
            ax1.axhline(cfl_maximum, color='darkslategray', linestyle=':', linewidth=1.2, label=f'Max: {cfl_maximum:.2f}', zorder=1)
            ax1.set_ylabel("CFL")
            ax1.set_title(f"CFL Over Time (target timestep: {self.dt})")
            ax1.legend(loc='best')
            ax1.grid(True)
            ax1.ticklabel_format(style='plain', axis='y', useOffset=False)

            # --- Plot Fixed Volume ---
            if has_fixed_volume:
                ax2.plot(fixed_volume_frames, linestyle='-', linewidth=2.5, color='red', label='Fixed Volume', zorder=2)
                ax2.axhline(fixed_volume_mean, color='darkgoldenrod', linestyle='--', linewidth=1.2, label=f'Mean: {fixed_volume_mean:.2f}', zorder=1)
                ax2.axhline(fixed_volume_median, color='sienna', linestyle='--', linewidth=1.2, label=f'Median: {fixed_volume_median:.2f}', zorder=1)
                ax2.axhline(fixed_volume_maximum, color='black', linestyle=':', linewidth=1.2, label=f'Min: {fixed_volume_minimum:.2f}', zorder=1)
                ax2.axhline(fixed_volume_maximum, color='black', linestyle=':', linewidth=1.2, label=f'Max: {fixed_volume_maximum:.2f}', zorder=1)
                ax2.set_ylabel("Fixed Volume")
                ax2.set_title("Fixed Volume Over Time")
                ax2.legend(loc='best')
                ax2.grid(True)
                    # Disable scientific notation and offset
                ax2.ticklabel_format(style='plain', axis='y', useOffset=False)

            # Common X label
            plt.xlabel("Frame")
            plt.tight_layout()
            plt.savefig(self.stats_dir / "CFL_and_Volume.png" if has_fixed_volume else "CFL.png", dpi=300)

    def finish(self):
        self.computeStats()
        
        if self.export_data:
            with open(self.stats_dir / "data.json", "w") as f:
                json.dump(self.data, f, indent=4)

        if self.export_videos:
            for index in self.tracked_grids_indeces:
                name = self.trackable_grids[index][0]


                frames_dir   = self.stats_dir / f"{name}_frames"
                pattern      = f"{name}_%04d.png"       # frame_0001.png, frame_0002.png …
                fps          = 30                       # playback frame‑rate
                crf          = 28                       # 0 (lossless) … 63 (worst). 28~30 ≈ YouTube HD
                bitrate      = "0"                      # Leave "0" for constrained‑quality mode
                output_file  = self.stats_dir / f"{name}.webm"

                cmd = [
                    "ffmpeg",
                    "-y",                               # overwrite existing output
                    "-framerate", str(fps),             # input fps
                    "-i", str(frames_dir / pattern),    # numbered‑file pattern
                    "-c:v", "libvpx-vp9",               # VP9 video codec
                    "-crf", str(crf),                   # quality target
                    "-b:v", bitrate,                    # needed even when "0"
                    "-pix_fmt", "yuv420p",              # 8‑bit 4:2:0 (widest support)
                    # Add an Opus audio track if desired (remove these two lines otherwise):
                    # "-f", "lavfi", "-i", "anullsrc=channel_layout=stereo:sample_rate=48000",
                    # "-c:a", "libopus", "-b:a", "128k",
                    output_file.as_posix(),
                ]

                try:
                    print("Running:", " ".join(cmd))
                    subprocess.run(cmd, check=True)
                except Exception as e:
                    print(f"could not export {name} to video becase {repr(e)}")

    def format_results(self, results, float_fmt="{:.6g}"):
        # Zeilen sammeln ──> [(name, metric, value_as_str), …]
        rows = []
        for name, stats in results.items():
            for metric, val in stats.items():
                val_str = float_fmt.format(val) if isinstance(val, (float, int)) else str(val)
                rows.append((name, metric, val_str))

        # Breite jeder Spalte ermitteln
        headers = ("Name", "Metric", "Value")
        col_w = [max(len(str(h)), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]

        # Hilfsfunktion für eine Zeile
        def make_row(items):
            return " | ".join(str(item).ljust(col_w[i]) for i, item in enumerate(items))

        # Tabelle zusammenbauen
        sep = "-+-".join("-" * w for w in col_w)
        lines = [make_row(headers), sep] + [make_row(r) for r in rows]
        return "\n".join(lines)

    def printStats(self):
        print(self.format_results(self.data["results"]))
