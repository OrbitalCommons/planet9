#!/usr/bin/env python3
"""Layer-2 coverage builder (design/03): alert store -> rubin_watch/lsst_coverage.json.

Derives LSST observed coverage from the alert parquet store. Each visit is
reconstructed by clustering detections on (band, exposure midpoint MJD); the
visit's footprint is painted as all HEALPix nside-64 pixels within the camera
radius of the detection centroid (`--mode dilated`, default) or as only the
pixels that actually contain detections (`--mode detections`, a strict lower
bound — appropriate once selection-B densities exist).

Depth bookkeeping: Light-SSO packets carry no per-visit limiting magnitude,
so `depth_single_median` records the labeled band fiducial for visited pixels
(design/07 note); the faintest detected magnitude per band is printed as a
sanity check that fiducials aren't overclaiming.

Usage:
  build_coverage.py [--store ~/.cache/p9-rubin-watch/alerts]
                    [--out rubin_watch/lsst_coverage.json]
                    [--mode dilated|detections]
"""

import argparse
import datetime
import json
import math
import pathlib
import sys

NSIDE = 64
NPIX = 12 * NSIDE * NSIDE
# LSST single-visit 5-sigma point-source fiducials (mag), labeled reference
# values (survey overview-paper scale); r matches the workspace's
# p9-2023-lsst-strategy SINGLE_VISIT_DEPTH_R.
BAND_FIDUCIAL = {"u": 23.9, "g": 24.8, "r": 24.3, "i": 23.9, "z": 23.3, "y": 22.1}
CAMERA_RADIUS_DEG = 1.75  # LSSTCam field of view is ~9.6 deg^2 (~3.5 deg wide)
VISIT_TIME_S = 2.0  # detections within this window + same band = one visit
MJD_UNIX_EPOCH = 40587.0


def now_utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def mjd_to_night(mjd: float) -> str:
    dt = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc) + datetime.timedelta(
        days=mjd - MJD_UNIX_EPOCH
    )
    return dt.strftime("%Y-%m-%d")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--store", type=pathlib.Path,
                    default=pathlib.Path.home() / ".cache" / "p9-rubin-watch" / "alerts")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("rubin_watch/lsst_coverage.json"))
    ap.add_argument("--mode", choices=["dilated", "detections"], default="dilated")
    args = ap.parse_args()

    import healpy as hp
    import numpy as np
    import pandas as pd
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    parts = sorted(args.store.rglob("part-*.parquet"))
    if not parts:
        print(f"no parquet under {args.store}", file=sys.stderr)
        return 2
    df = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    df = df[df.survey == "lsst"].drop_duplicates(subset="alert_id")
    if df.empty:
        print("no LSST records in store", file=sys.stderr)
        return 2

    # --- reconstruct visits: cluster on (band, quantized mjd) -------------
    df = df.sort_values("mjd")
    df["visit_key"] = df.band + "|" + (df.mjd * 86400 / VISIT_TIME_S).round().astype("int64").map(str)
    visits = df.groupby("visit_key").agg(
        band=("band", "first"), mjd=("mjd", "median"),
        ra=("ra", "median"), dec=("dec", "median"), n_det=("alert_id", "size"),
    )

    bands = sorted(set(visits.band) & set(BAND_FIDUCIAL))
    n_visits = {b: np.zeros(NPIX, dtype=np.int32) for b in bands}
    last_mjd = {b: np.zeros(NPIX, dtype=np.float64) for b in bands}

    for _, v in visits.iterrows():
        if v.band not in n_visits:
            continue
        if args.mode == "dilated":
            vec = hp.ang2vec(v.ra, v.dec, lonlat=True)
            pix = hp.query_disc(NSIDE, vec, math.radians(CAMERA_RADIUS_DEG),
                                inclusive=True)
        else:
            sel = df[df.visit_key == v.name]
            pix = np.unique(hp.ang2pix(NSIDE, sel.ra.to_numpy(),
                                       sel.dec.to_numpy(), lonlat=True))
        n_visits[v.band][pix] += 1
        last_mjd[v.band][pix] = np.maximum(last_mjd[v.band][pix], v.mjd)

    # --- flags -------------------------------------------------------------
    lon, lat = hp.pix2ang(NSIDE, np.arange(NPIX), lonlat=True)
    gal_b = SkyCoord(ra=lon * u.deg, dec=lat * u.deg, frame="icrs").galactic.b.deg
    crowding = np.abs(gal_b) < 10.0
    template_epoch = np.zeros(NPIX, dtype=bool)  # year-1: no survey-data templates yet

    # --- assemble ----------------------------------------------------------
    out = {
        "schema_version": 1,
        "generated_utc": now_utc(),
        "healpix": {"nside": NSIDE, "ordering": "ring"},
        "window": {
            "first_night": mjd_to_night(float(df.mjd.min())),
            "last_night": mjd_to_night(float(df.mjd.max())),
        },
        "source": f"fink_alert_metadata ({args.mode}; camera_radius_deg="
                  f"{CAMERA_RADIUS_DEG if args.mode == 'dilated' else None})",
        "n_input_records": int(len(df)),
        "n_reconstructed_visits": int(len(visits)),
        "band_fiducial_depth": {b: BAND_FIDUCIAL[b] for b in bands},
        # Sparse per-band encoding: aligned lists over covered pixels only
        # (single-visit depth is the band fiducial for every covered pixel,
        # so it is carried once above, not per pixel).
        "bands": {
            b: {
                "pixels": np.nonzero(n_visits[b])[0].tolist(),
                "n_visits": n_visits[b][n_visits[b] > 0].tolist(),
                "last_visit_mjd": [
                    round(float(x), 3) for x in last_mjd[b][n_visits[b] > 0]
                ],
            }
            for b in bands
        },
        "flags": {
            "template_epoch_pixels": np.nonzero(template_epoch)[0].tolist(),
            "crowding_pixels": np.nonzero(crowding)[0].tolist(),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, separators=(",", ":")) + "\n")

    covered = np.zeros(NPIX, dtype=bool)
    for b in bands:
        covered |= n_visits[b] > 0
    pix_area = 4 * math.pi / NPIX * (180 / math.pi) ** 2
    print(
        f"wrote {args.out} ({args.out.stat().st_size/1e6:.2f} MB): "
        f"{len(visits)} visits, bands {bands}, "
        f"{int(covered.sum())} pixels covered ({covered.sum()*pix_area:.0f} deg2, "
        f"{100*covered.mean():.2f}% of sky)"
    )
    for b in bands:
        sel = df[df.band == b]
        print(f"  {b}: {int((n_visits[b]>0).sum())} pix, faintest detection "
              f"{sel.psf_mag.max():.2f} vs fiducial {BAND_FIDUCIAL[b]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
