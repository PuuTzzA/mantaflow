import os
from manta import *
from pathlib import Path
import subprocess
import shutil
import json

class Data_collectior:
    def __init__(self, title="no_title_specified", base_dir="../analysis/experiments/", params=None, export_data=True, export_images=False):
        self.title = title
        self.base_dir = Path(base_dir).expanduser().resolve()
        self.stats_dir = self.base_dir / f"{self.title}_stats"

        if params is not None:
            self.dim = params["dimension"]
            self.res = [params["resolutionX"], params["resolutionY"], params["resolutionZ"]]
            self.cfl = params["CFL"]
            self.max_time = params["max_time"]
            self.fps = params["fps"]
        
        self.export_data = export_data
        self.export_images = export_images

    def init(self):
        shutil.rmtree(self.stats_dir, ignore_errors=True)
        self.stats_dir.mkdir(parents=True, exist_ok=True)
        
        self.data = {}
        self.data["title"] = self.title
        self.data["dim"] = self.dim
        self.data["res"] = self.res
        self.data["cfl"] = self.cfl
        self.data["max_time"] = self.max_time
        self.data["fps"] = self.fps
        self.data["frame_data"] = {}

        self.tracked_grids = set()
        
        if (self.export_images):
           (self.stats_dir / "frames").mkdir(parents=True, exist_ok=True)

    def step(self, solver, tracked_grids, flags, gui=None, windowSize=[800, 800], camPos=[0, 0, -1.3]):
        self.data["frame_data"][str(solver.frame).zfill(4)] = {}

        for [grid, name] in tracked_grids:
            self.data["frame_data"][str(solver.frame).zfill(4)][name] = json.loads(realGridStats(grid=grid, flags=flags))
            self.tracked_grids.add(name)

        if self.export_images and gui is not None:
            gui.windowSize(windowSize[0], windowSize[1])
            gui.setCamPos(camPos[0], camPos[1], camPos[2])
            gui.screenshot(str(self.stats_dir / "frames" / f"frame_{str(solver.frame).zfill(4)}.png"))

    def finish(self):
        if self.export_data:
            with open(self.stats_dir / "data.json", "w") as f:
                json.dump(self.data, f, indent=4)

        if self.export_images:
            frames_dir   = self.stats_dir / "frames"
            pattern      = "frame_%04d.png"        # frame_0001.png, frame_0002.png …
            fps          = self.fps                # playback frame‑rate
            crf          = 28                      # 0 (lossless) … 63 (worst). 28~30 ≈ YouTube HD
            bitrate      = "0"                     # Leave "0" for constrained‑quality mode
            output_file  = self.stats_dir / "simulation.webm"

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

            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)

    def printStats(self):
        frames = self.data["frame_data"]
        if not frames:
            print("No frame data collected")
            return

        first_frame_key = min(frames.keys())
        name_width = max(len(name) for name in self.tracked_grids) if self.tracked_grids else 0

        for grid_name in sorted(self.tracked_grids):
            start_sum = self.data["frame_data"][first_frame_key][grid_name]["sum"]
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
                # starting value = value in the earliest frame (ideally "0000")
                if f_key == first_frame_key:
                    start_sum = s
                # track min / max over all frames
                min_sum = min(min_sum, s)
                max_sum = max(max_sum, s)

            if not found_any:
                print(f"[{grid_name:<{name_width}}] sum: N/A (no values found)")
                continue

            # If the earliest frame lacked this grid, start_sum may still be None
            start_display = f"{start_sum:,.2f}" if start_sum is not None else "N/A"
            print(f"[{grid_name:<{name_width}}] sum: start = {start_display}, "
                  f"min = {min_sum:,.2f}, max = {max_sum:,.2f}")