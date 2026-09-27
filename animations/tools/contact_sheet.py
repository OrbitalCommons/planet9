#!/usr/bin/env python3
"""Render one scene at low quality and tile evenly spaced frames into a PNG,
so a scene's layout can be reviewed at a glance.

Usage (from animations/):
    python tools/contact_sheet.py scenes/papers/p9_2021_ztf/scene.py Ztf2021 [--frames 12] [--out sheet.png]

Frames are sampled at the midpoint of equal time slices, labelled with their
timestamp, and tiled 3 across. Needs manim and ffmpeg on PATH.
"""
import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("scene")
    ap.add_argument("--frames", type=int, default=12)
    ap.add_argument("--out", default=None)
    ap.add_argument("--quality", default="l", choices=["l", "m", "h"])
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix="p9sheet_")
    env = {**os.environ, "PYTHONPATH": HERE,
           "PATH": os.path.expanduser("~/.local/bin") + ":" + os.environ.get("PATH", "")}
    r = subprocess.run(["manim", "render", f"-q{args.quality}", "--media_dir", work,
                        os.path.join(HERE, args.file), args.scene],
                       cwd=HERE, env=env, capture_output=True, text=True)
    mp4 = None
    for root, _, files in os.walk(os.path.join(work, "videos")):
        if "partial_movie_files" in root:
            continue
        for f in files:
            if f == args.scene + ".mp4":
                mp4 = os.path.join(root, f)
    if r.returncode != 0 or mp4 is None:
        sys.stderr.write(r.stdout[-3000:] + "\n" + r.stderr[-3000:] + "\n")
        sys.exit(f"render failed: {args.scene}")

    dur = float(subprocess.run(
        ["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0", mp4],
        capture_output=True, text=True, check=True).stdout.strip())
    n = args.frames
    shots = []
    for k in range(n):
        t = dur * (k + 0.5) / n
        png = os.path.join(work, f"f{k:02d}.png")
        label = f"{t:5.1f}s".replace(":", "\\:")
        subprocess.run(
            ["ffmpeg", "-y", "-v", "error", "-ss", f"{t:.3f}", "-i", mp4, "-frames:v", "1",
             "-vf", f"scale=640:-1,drawtext=text='{label}':x=8:y=8:fontsize=18:"
                    "fontcolor=white:box=1:boxcolor=black@0.6",
             png], check=True)
        shots.append(png)
    out = args.out or os.path.join(work, f"{args.scene}_sheet.png")
    cols = 3
    rows = (n + cols - 1) // cols
    inputs = []
    for s in shots:
        inputs += ["-i", s]
    layout = "|".join(f"{(k % cols) * 640}_{(k // cols) * 360}" for k in range(n))
    subprocess.run(["ffmpeg", "-y", "-v", "error", *inputs, "-filter_complex",
                    f"xstack=inputs={n}:layout={layout}:fill=black", "-frames:v", "1", out],
                   check=True)
    print(f"{out}  ({dur:.1f}s, {n} frames, {cols}x{rows})")


if __name__ == "__main__":
    main()
