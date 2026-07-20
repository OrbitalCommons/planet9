#!/usr/bin/env python3
"""Export linker inputs (design/05) from the alert store.

Per HEALPix nside-8 tile with any detections in the window, writes a
LinkerInput JSON: the tile's detections (unassociated only, unless
--include-known for smoke tests on selection-A data) plus the reconstructed
visit list (same band+epoch clustering as the coverage builder) — the visit
list drives the injection/recovery calibration.

Usage:
  export_linker_input.py [--store DIR] [--outdir DIR]
                         [--mjd-start X --mjd-end Y] [--include-known]
"""

import argparse
import json
import math
import pathlib
import sys

TILE_NSIDE = 8
CAMERA_RADIUS_DEG = 1.75
VISIT_TIME_S = 2.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--store", type=pathlib.Path,
                    default=pathlib.Path.home() / ".cache" / "p9-rubin-watch" / "alerts")
    ap.add_argument("--outdir", type=pathlib.Path,
                    default=pathlib.Path.home() / ".cache" / "p9-rubin-watch" / "linker")
    ap.add_argument("--mjd-start", type=float, default=None)
    ap.add_argument("--mjd-end", type=float, default=None)
    ap.add_argument("--include-known", action="store_true",
                    help="keep ss-associated detections (selection-A smoke tests)")
    args = ap.parse_args()

    import healpy as hp
    import numpy as np
    import pandas as pd

    parts = sorted(args.store.rglob("part-*.parquet"))
    if not parts:
        print(f"no parquet under {args.store}", file=sys.stderr)
        return 2
    df = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    df = df[df.survey == "lsst"].drop_duplicates(subset="alert_id")
    if args.mjd_start is not None:
        df = df[df.mjd >= args.mjd_start]
    if args.mjd_end is not None:
        df = df[df.mjd <= args.mjd_end]
    if df.empty:
        print("no records in window", file=sys.stderr)
        return 2

    # Visit reconstruction over the FULL window (visits cross tile edges).
    df = df.sort_values("mjd")
    df["visit_key"] = df.band + "|" + (df.mjd * 86400 / VISIT_TIME_S).round().astype("int64").map(str)
    visits = df.groupby("visit_key").agg(
        band=("band", "first"), mjd=("mjd", "median"),
        ra=("ra", "median"), dec=("dec", "median"),
    ).reset_index(drop=True)

    keep = df if args.include_known else df[df.ss_name.isna()]
    if keep.empty:
        print("no unassociated detections in window (selection-B not pulled yet?); "
              "use --include-known for smoke tests", file=sys.stderr)
        return 2
    keep = keep.copy()
    keep["tile"] = hp.ang2pix(TILE_NSIDE, keep.ra.to_numpy(), keep.dec.to_numpy(),
                              lonlat=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    written = 0
    for tile, sel in keep.groupby("tile"):
        # Visits whose cone could touch the tile: center within tile radius
        # (~4 deg for nside 8) + camera radius of any tile detection centroid.
        c_ra, c_dec = float(sel.ra.median()), float(sel.dec.median())
        sep = hp.rotator.angdist(
            np.vstack([visits.ra.to_numpy(), visits.dec.to_numpy()]),
            np.array([c_ra, c_dec]), lonlat=True,
        )
        vsel = visits[sep < math.radians(6.0 + CAMERA_RADIUS_DEG)]
        out = {
            "tile": int(tile),
            "mjd_start": float(keep.mjd.min()),
            "mjd_end": float(keep.mjd.max()),
            "detections": [
                {
                    "alert_id": int(r.alert_id),
                    "mjd": float(r.mjd),
                    "ra": float(r.ra),
                    "dec": float(r.dec),
                    "band": str(r.band),
                    "psf_mag": None if pd.isna(r.psf_mag) else float(r.psf_mag),
                    "sigma_arcsec": None,
                }
                for r in sel.itertuples()
            ],
            "visits": [
                {
                    "mjd": float(v.mjd),
                    "band": str(v.band),
                    "ra": float(v.ra),
                    "dec": float(v.dec),
                    "radius_deg": CAMERA_RADIUS_DEG,
                }
                for v in vsel.itertuples()
            ],
        }
        path = args.outdir / f"tile-{int(tile):04}.json"
        path.write_text(json.dumps(out) + "\n")
        written += 1
        print(f"  tile {int(tile)}: {len(sel)} detections, {len(vsel)} visits -> {path.name}")
    print(f"wrote {written} tile inputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
