import os
from manta import *
from pathlib import Path
import subprocess
import shutil
import json

class Data_collectior:
    def __init__(self, title="no_title_specified", base_dir="../analysis/experiments/", params=None, export_data=True, export_images=False, trackable_grid_names=[], tracked_grids_indeces=[]):
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
        self.trackable_grid_names=trackable_grid_names
        self.tracked_grids_indeces = tracked_grids_indeces

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

        self.tracked_grids = set()
        self.last_frame = -1
        
        if (self.export_images):
           for index in self.tracked_grids_indeces:
            (self.stats_dir / f"{self.trackable_grid_names[index]}_frames").mkdir(parents=True, exist_ok=True)

    def step(self, solver, tracked_grids, flags, vel, gui=None, windowSize=[800, 800], camPos=[0, 0, -1.3]):
        if (self.last_frame == solver.frame):
            return

        self.last_frame = solver.frame
        self.data["frame_data"][str(solver.frame).zfill(4)] = {}
        self.data["frame_data"][str(solver.frame).zfill(4)]["cfl"] = vel.getMaxAbs() * solver.timestep

        for [grid, name] in tracked_grids:
            self.data["frame_data"][str(solver.frame).zfill(4)][name] = json.loads(realGridStats(grid=grid, flags=flags))
            self.tracked_grids.add(name)

        if self.export_images and gui is not None:
            gui.windowSize(windowSize[0], windowSize[1])
            gui.setCamPos(camPos[0], camPos[1], camPos[2])

            for i in range(len(self.trackable_grid_names)):
                if i in self.tracked_grids_indeces:
                    name = self.trackable_grid_names[i]
                    gui.screenshot(str(self.stats_dir / f"{name}_frames" / f"{name}_{str(solver.frame).zfill(4)}.png"))
                
                gui.nextRealGrid()
                gui.update()

    def computeStats(self):
        frames = self.data["frame_data"]
        if not frames:
            print("No frame data collected")
            return

        self.data["results"] = {}

        min_cfl = float("inf")
        max_cfl = -float("inf")
        for f_key, grids in frames.items():
            try:
                current_cfl = grids["cfl"]
            except (KeyError):
                continue
            min_cfl = min(min_cfl, current_cfl)
            max_cfl = max(max_cfl, current_cfl)
        
        self.data["results"]["cfl"] = {}
        self.data["results"]["cfl"]["minCfl"] = min_cfl
        self.data["results"]["cfl"]["maxCfl"] = max_cfl

        for grid_name in sorted(self.tracked_grids):
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
            with open(self.stats_dir / "data.json", "w") as f:
                json.dump(self.data, f, indent=4)

        if self.export_images:
            for index in self.tracked_grids_indeces:
                name = self.trackable_grid_names[index]


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
        """
        results  : dict wie data["results"]
        float_fmt: Format‑String für Gleitkommazahlen (default 6 signifikante Stellen)

        Gibt einen formatierten String zurück, der sich einfach ausdrucken lässt.
        """
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
