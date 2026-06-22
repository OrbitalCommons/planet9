#!/usr/bin/env python3
"""Render every scene in manifest.yaml (in parallel) and concatenate into one film.

Usage:
    PYTHONPATH=. python assembler.py [--quality l|m|h] [--jobs N] [--render-only]

Each scene renders into its own isolated --media_dir (so parallel manim processes
never collide on partial-movie / TeX caches), then the per-scene mp4s are
concatenated in manifest order into animations/out/planet9_film.mp4.
Needs manim on PATH plus ffmpeg.
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p9_manim import content  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(HERE, "out")
WORK = "/tmp/p9film"
QMAP = {"l": "480p15", "m": "720p30", "h": "1080p60", "k": "2160p60"}
ENV = {**os.environ, "PYTHONPATH": HERE,
       "PATH": os.path.expanduser("~/.local/bin") + ":" + os.environ.get("PATH", "")}


def manifest_entries():
    with open(os.path.join(HERE, "manifest.yaml")) as fh:
        m = yaml.safe_load(fh)
    seen, ordered = set(), []
    for group in ("preface", "papers"):
        for e in m.get(group) or []:
            key = (e["file"], e["scene"])
            if key not in seen:
                seen.add(key)
                ordered.append(e)
    return ordered


def with_companions(entries):
    """Interleave each scene's auto-generated discussion/learning beat after it."""
    out = []
    for e in entries:
        out.append(e)
        f = e["file"]
        if "/preface/" in f and e["scene"] in content.PREFACE_LEARNING:
            out.append({"file": "scenes/discussion.py", "scene": "Learn_" + e["scene"]})
        elif "/papers/" in f:
            crate = os.path.basename(os.path.dirname(f)).replace("_", "-")
            if crate in content.PAPER_TEXT:
                out.append({"file": "scenes/discussion.py",
                            "scene": "Discuss_" + crate.replace("-", "_")})
    return out


def scene_mp4(idx, entry, qdir):
    stem = os.path.splitext(os.path.basename(entry["file"]))[0]
    return os.path.join(WORK, str(idx), "videos", stem, qdir, entry["scene"] + ".mp4")


def render(idx, entry, quality, qdir):
    target = scene_mp4(idx, entry, qdir)
    if os.path.exists(target):
        return idx, entry, True, ""
    media_dir = os.path.join(WORK, str(idx))
    r = subprocess.run(
        ["manim", "render", f"-q{quality}", "--media_dir", media_dir,
         os.path.join(HERE, entry["file"]), entry["scene"]],
        cwd=HERE, env=ENV, capture_output=True, text=True,
    )
    ok = r.returncode == 0 and os.path.exists(target)
    return idx, entry, ok, (r.stderr[-400:] if not ok else "")


def concat(mp4s, out_path):
    listfile = os.path.join(OUT, "_concat.txt")
    with open(listfile, "w") as fh:
        for p in mp4s:
            fh.write(f"file '{p}'\n")
    subprocess.run(["ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", listfile,
                    "-c", "copy", out_path], check=True, capture_output=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quality", default="l", choices=list(QMAP))
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    qdir = QMAP[args.quality]
    entries = with_companions(manifest_entries())
    print(f"{len(entries)} scenes (incl. discussion/learning), {args.jobs} parallel jobs", flush=True)

    results = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(render, i, e, args.quality, qdir) for i, e in enumerate(entries)]
        done = 0
        for fut in as_completed(futs):
            idx, entry, ok, err = fut.result()
            results[idx] = (entry, ok)
            done += 1
            status = "ok" if ok else "FAIL"
            print(f"[{done}/{len(entries)}] {status}: {entry['scene']}", flush=True)
            if not ok:
                print(f"    {err}", file=sys.stderr, flush=True)

    mp4s, failed = [], []
    for idx in range(len(entries)):
        entry, ok = results[idx]
        if ok:
            mp4s.append(scene_mp4(idx, entry, qdir))
        else:
            failed.append(entry["scene"])
    print(f"\nrendered {len(mp4s)}/{len(entries)}; failed: {failed}", flush=True)

    if args.render_only:
        return
    out_path = os.path.join(OUT, "planet9_film.mp4")
    concat(mp4s, out_path)
    print(f"FILM: {out_path}  ({len(mp4s)} scenes)", flush=True)


if __name__ == "__main__":
    main()
