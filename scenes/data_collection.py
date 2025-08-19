import os
from manta import *
from pathlib import Path
import subprocess
import shutil
import json
import numpy as np
import matplotlib.pyplot as plt
import math
from graph_creation import *

class Data_collectior:
    def __init__(self, title="no_title_specified", base_dir="../exports/experiments/", params=None, export_data=True, export_images=False, export_videos=False, export_vdbs=False, trackable_grid_names=[], tracked_grids_indeces=[], image_grids_indeces=[], graph_grids=[]):
        self.title = title
        self.base_dir = Path(base_dir).expanduser().resolve() / self.title

        if params is not None:
            self.dim = params["dimension"]
            self.res = [params["resolutionX"], params["resolutionY"], params["resolutionZ"]]
            self.maxCfl = params["maxCFL"]
            self.max_time = params["max_time"]
            self.dt = params["dt"]
        
        self.export_data = export_data
        self.export_images = export_images
        self.export_videos = export_videos
        self.export_vdbs = export_vdbs
        self.trackable_grids=trackable_grid_names
        self.tracked_grids_indeces = tracked_grids_indeces
        self.image_grids_indeces = image_grids_indeces
        self.graph_grids = graph_grids

    def init(self):
        #if self.export_data or self.export_images:
        #    shutil.rmtree(self.base_dir, ignore_errors=True)
        if self.export_data or self.export_images or self.export_videos or self.export_vdbs:
            self.base_dir.mkdir(parents=True, exist_ok=True)
        
        self.data = {}
        self.data["title"] = self.title
        self.data["dim"] = self.dim
        self.data["res"] = self.res
        self.data["maxCfl"] = self.maxCfl
        self.data["max_time"] = self.max_time
        self.data["dt"] = self.dt
        self.data["frame_data"] = {}

        self.current_frame = 0
        self.last_framerate_frame = -1
        
        if (self.export_images):
           for index in self.image_grids_indeces:
            (self.base_dir / f"{self.trackable_grids[index][0]}_frames").mkdir(parents=True, exist_ok=True)
        
        if self.export_vdbs:
            (self.base_dir / "vdbs").mkdir(parents=True, exist_ok=True)

    def step(self, solver, flags, maxVel, gui=None, windowSize=[1000, 1000], camPos=[0, 0, -1.2], objects=[]):
        #self.current_frame = math.floor(solver.timeTotal)
        self.data["frame_data"][str(self.current_frame).zfill(4)] = {}
        self.data["frame_data"][str(self.current_frame).zfill(4)]["cfl"] = maxVel * solver.timestep
        self.data["frame_data"][str(self.current_frame).zfill(4)]["dt"] = solver.timestep

        for i in range(len(self.trackable_grids)):
            if i in self.tracked_grids_indeces:    
                name = self.trackable_grids[i][0]              
                grid = self.trackable_grids[i][1]
                self.data["frame_data"][str(self.current_frame).zfill(4)][name] = json.loads(realGridStats(grid=grid, flags=flags))

        if self.export_images and gui is not None and (int(solver.timeTotal * 0.4) != self.last_framerate_frame): #and (self.last_framerate_frame != solver.frame):
            gui.windowSize(windowSize[0], windowSize[1])
            gui.setCamPos(camPos[0], camPos[1], camPos[2])
            gui.update()
            for i in range(len(self.trackable_grids)):
                if i in self.image_grids_indeces:
                    name = self.trackable_grids[i][0]
                    #gui.screenshot(str(self.base_dir / f"{name}_frames" / f"{name}_{str(math.floor(solver.frame)).zfill(4)}.png"))
                    gui.screenshot(str(self.base_dir / f"{name}_frames" / f"{name}_{str(math.floor(int(solver.timeTotal * 0.4))).zfill(4)}.png"))
                    #gui.screenshot(str(self.base_dir / f"{name}_frames" / f"{name}_{str(math.floor(self.current_frame)).zfill(4)}.png"))
                    
                gui.nextRealGrid()
                gui.update()
        
        if self.export_vdbs and (self.last_framerate_frame != solver.frame):
            # note: when saving pdata fields, they must be accompanied by and listed before their parent pp
            save( name=str(self.base_dir / "vdbs" / f"{self.title}_{solver.frame}.vdb"), objects=objects )

        self.current_frame += 1
        self.last_framerate_frame = solver.frame
        self.last_framerate_frame = int(solver.timeTotal * 0.4)

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

            if not found_any:
                continue

            self.data["results"][grid_name] = {}
            self.data["results"][grid_name]["minSum"] = min_sum
            self.data["results"][grid_name]["maxSum"] = max_sum

    def finish(self):
        self.computeStats()
        
        if self.export_data:
            with open(self.base_dir / "data.json", "w") as f:
                json.dump(self.data, f, indent=4)

            create_combined_graph(data_array=[self.data], data_names=[""], interested_fields=self.graph_grids, 
                                  title=self.title, include_cfl=True, include_dt=True, include_extra_stats=True, export_path=self.base_dir / "graph.png")

        if self.export_videos:
            for index in self.image_grids_indeces:
                name = self.trackable_grids[index][0]

                frames_dir   = self.base_dir / f"{name}_frames"
                pattern      = f"{name}_%04d.png"       # frame_0001.png, frame_0002.png …
                fps          = 30                       # playback frame‑rate
                crf          = 28                       # 0 (lossless) … 63 (worst). 28~30 ≈ YouTube HD
                bitrate      = "0"                      # Leave "0" for constrained‑quality mode
                output_file  = self.base_dir / f"{name}.webm"

                # Title overlay (drawtext filter)
                drawtext_filter = f"drawtext=text='{self.title}':fontsize=24:fontcolor=white:x=10:y=10:borderw=2"

                cmd = [
                    "ffmpeg",
                    "-y",                               # overwrite existing output
                    "-framerate", str(fps),             # input fps
                    "-i", str(frames_dir / pattern),    # numbered‑file pattern
                    "-vf", drawtext_filter,             # title overlay
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
