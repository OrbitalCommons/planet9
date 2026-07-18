#!/usr/bin/env python3
"""Fink Data Transfer pull for the Rubin watch (Layer 3, design/04).

Registration landed 2026-07-17 (see rubin_watch/RUNBOOK-fink.md for the
verified environment: uv venv at ~/.venvs/fink, fink-client >= 11,
servers kafka-{lsst,ztf}.fink-broker.org:24499). Credentials live in
~/.finkclient/{lsst,ztf}_credentials.yml; this script refuses to run
without them rather than half-working.

Usage:
  fink_pull.py --check-creds
  fink_pull.py --topic ftransfer_lsst_... --outdir ~/.cache/p9-rubin-watch/alerts
  fink_pull.py --prune --outdir ... [--retention-days 120]

The Data Transfer job itself is created in the portal
(https://lsst.fink-portal.org/download); this script consumes the resulting
topic to parquet partitions night=YYYYMMDD/tile=NNN (schema: design/07) and
enforces the volume alarm (design/08).
"""

import argparse
import pathlib
import shutil
import sys

CRED_PATH = pathlib.Path.home() / ".finkclient" / "lsst_credentials.yml"
# design/04 volume budget: worst case 1e5 selection-B records/night; the ops
# alarm trips at 3x.
NIGHTLY_RECORD_BUDGET = 100_000
ALARM_FACTOR = 3


def check_creds() -> bool:
    if CRED_PATH.exists():
        print(f"credentials present: {CRED_PATH}")
        return True
    print(
        f"no Fink LSST credentials at {CRED_PATH}\n"
        "  1. fill the subscription form (fink-client README) with the LSST streams\n"
        "  2. store username/group in 1Password; op inject into env\n"
        "  3. pip install fink-client && fink_client_register -survey lsst \\\n"
        "       -username $FINK_USERNAME -group_id $FINK_GROUP -servers <from email>",
        file=sys.stderr,
    )
    return False


def prune(outdir: pathlib.Path, retention_days: int) -> None:
    """Drop night= partitions older than the retention window."""
    nights = sorted(p for p in outdir.glob("night=*") if p.is_dir())
    if len(nights) <= retention_days:
        print(f"{len(nights)} night partitions <= retention {retention_days}; nothing to prune")
        return
    for p in nights[: len(nights) - retention_days]:
        print(f"pruning {p}")
        shutil.rmtree(p)


# ---------------------------------------------------------------- extraction

MJD_UNIX_EPOCH = 40587.0  # MJD of 1970-01-01


def _mjd_to_night(mjd: float) -> str:
    """UTC night label YYYYMMDD for a mid-exposure MJD."""
    import datetime

    dt = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc) + datetime.timedelta(
        days=mjd - MJD_UNIX_EPOCH
    )
    return dt.strftime("%Y%m%d")


def extract_record(alert: dict, survey: str, topic: str) -> dict | None:
    """Map one alert packet to the design/07 parquet record (columns doc'd
    there). ZTF and LSST packets differ; both mappings live here so the ZTF
    SSO testbed exercises the same code path the LSST transfers will use.
    Returns None when the packet lacks the positional minimum."""
    if survey == "ztf":
        c = alert.get("candidate") or {}
        if c.get("ra") is None or c.get("jd") is None:
            return None
        return {
            "alert_id": int(c.get("candid") or 0),
            "dia_source_id": int(c.get("candid") or 0),
            "dia_object_id": None,  # ZTF has no numeric diaObjectId; objectId below
            "ss_object_id": None if c.get("ssnamenr") in (None, "null") else 1,
            "ss_name": None if c.get("ssnamenr") in (None, "null") else str(c.get("ssnamenr")),
            "object_id": str(alert.get("objectId")),
            "mjd": float(c["jd"]) - 2400000.5,
            "ra": float(c["ra"]),
            "dec": float(c["dec"]),
            "ra_err": None,
            "dec_err": None,
            "band": str(c.get("fid")),
            "psf_mag": None if c.get("magpsf") is None else float(c["magpsf"]),
            "psf_mag_err": None if c.get("sigmapsf") is None else float(c["sigmapsf"]),
            "reliability": None if c.get("drb") is None else float(c["drb"]),
            "visit": None,
            "detector": None,
            "fink_class": alert.get("fink_class"),
            "topic": topic,
            "survey": survey,
        }
    # LSST packet (diaSource-rooted; field names per the Rubin alert schema
    # as served by Fink — verified against the first Data Transfer pull).
    d = alert.get("diaSource") or {}
    if d.get("ra") is None or d.get("midpointMjdTai") is None:
        return None
    return {
        "alert_id": int(alert.get("alertId") or d.get("diaSourceId") or 0),
        "dia_source_id": int(d.get("diaSourceId") or 0),
        "dia_object_id": None if d.get("diaObjectId") is None else int(d["diaObjectId"]),
        "ss_object_id": None if d.get("ssObjectId") in (None, 0) else int(d["ssObjectId"]),
        "ss_name": None,
        "object_id": None if d.get("diaObjectId") is None else str(d["diaObjectId"]),
        "mjd": float(d["midpointMjdTai"]),
        "ra": float(d["ra"]),
        "dec": float(d["dec"]),
        "ra_err": None if d.get("raErr") is None else float(d["raErr"]),
        "dec_err": None if d.get("decErr") is None else float(d["decErr"]),
        "band": str(d.get("band")),
        "psf_mag": None if d.get("psfFlux") is None else _njy_to_mag(float(d["psfFlux"])),
        "psf_mag_err": None,
        "reliability": None if d.get("reliability") is None else float(d["reliability"]),
        "visit": None if d.get("visit") is None else int(d["visit"]),
        "detector": None if d.get("detector") is None else int(d["detector"]),
        "fink_class": alert.get("fink_class"),
        "topic": topic,
        "survey": survey,
    }


def _njy_to_mag(flux_njy: float) -> float | None:
    """AB magnitude from nJy PSF flux (LSST alert fluxes are nJy)."""
    if flux_njy <= 0:
        return None
    return -2.5 * __import__("math").log10(flux_njy * 1e-9) + 8.90


def write_partitioned(records: list[dict], outdir: pathlib.Path) -> dict:
    """Write records to night=YYYYMMDD/tile=NNN parquet partitions (nside-8
    ring healpix). Returns per-night counts; enforces the volume alarm."""
    import astropy.units as u
    import pandas as pd
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=8, order="ring")
    df = pd.DataFrame.from_records(records)
    df["night"] = df["mjd"].map(_mjd_to_night)
    df["tile"] = hp.lonlat_to_healpix(
        df["ra"].to_numpy() * u.deg, df["dec"].to_numpy() * u.deg
    ).astype("int32")

    counts = df.groupby("night").size().to_dict()
    for night, n in counts.items():
        if n > NIGHTLY_RECORD_BUDGET * ALARM_FACTOR:
            raise SystemExit(
                f"VOLUME ALARM: night {night} has {n} records "
                f"(> {ALARM_FACTOR}x budget {NIGHTLY_RECORD_BUDGET}); aborting before write"
            )
    for (night, tile), part in df.groupby(["night", "tile"]):
        pdir = outdir / f"night={night}" / f"tile={tile}"
        pdir.mkdir(parents=True, exist_ok=True)
        # Idempotent-ish append: one file per ingest batch, deduped on read
        # via alert_id (design/04 replay note).
        n_existing = len(list(pdir.glob("part-*.parquet")))
        part.drop(columns=["night", "tile"]).to_parquet(
            pdir / f"part-{n_existing:04}.parquet", index=False
        )
    return counts


def poll_livestream(survey: str, limit: int, outdir: pathlib.Path, timeout: float) -> int:
    """Consume up to `limit` alerts from the registered livestream topics and
    land them in the partitioned store. The ZTF SSO topics are the testbed;
    the same path serves LSST when its SSO streams appear."""
    if not check_creds():
        return 2
    from fink_client.configuration import load_credentials
    from fink_client.consumer import AlertConsumer

    creds = load_credentials(survey=survey)
    config = {
        "username": creds["username"],
        "bootstrap.servers": creds["servers"],
        "group.id": creds["group_id"],
    }
    consumer = AlertConsumer(creds["mytopics"], config, survey=survey)
    try:
        out = consumer.consume(num_alerts=limit, timeout=timeout)
    finally:
        consumer.close()
    records = []
    for topic, alert, _key in out:
        if topic is None or alert is None:
            continue
        r = extract_record(alert, survey, topic)
        if r is not None:
            records.append(r)
    if not records:
        print("no alerts within timeout (queue empty?)", file=sys.stderr)
        return 0
    counts = write_partitioned(records, outdir)
    print(f"landed {len(records)} records: {counts} -> {outdir}")
    return 0


def pull(topic: str, outdir: pathlib.Path) -> int:
    """Consume a Data Transfer topic (created in the portal) into the store."""
    if not check_creds():
        return 2
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        # fink_datatransfer is the reference consumer for ftransfer topics
        # (handles the sidecar _schema topic); we re-partition its output.
        rc = subprocess.run(
            [
                str(pathlib.Path.home() / ".venvs/fink/bin/fink_datatransfer"),
                "-survey",
                "lsst",
                "-topic",
                topic,
                "-outdir",
                tmp,
                "-nconsumers",
                "4",
                "--verbose",
            ],
            check=False,
        ).returncode
        if rc != 0:
            print(f"fink_datatransfer exited {rc}", file=sys.stderr)
            return 2
        import pandas as pd

        frames = [pd.read_parquet(p) for p in pathlib.Path(tmp).rglob("*.parquet")]
        if not frames:
            print("transfer produced no parquet", file=sys.stderr)
            return 2
        raw = pd.concat(frames, ignore_index=True)
        records = [
            r
            for r in (
                extract_record(row, "lsst", topic) for row in raw.to_dict("records")
            )
            if r is not None
        ]
        counts = write_partitioned(records, outdir)
        print(f"landed {len(records)} of {len(raw)} rows: {counts} -> {outdir}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-creds", action="store_true")
    ap.add_argument("--topic")
    ap.add_argument("--outdir", type=pathlib.Path,
                    default=pathlib.Path.home() / ".cache" / "p9-rubin-watch" / "alerts")
    ap.add_argument("--prune", action="store_true")
    ap.add_argument("--retention-days", type=int, default=120)
    ap.add_argument("--poll-livestream", action="store_true",
                    help="consume registered livestream topics into the store")
    ap.add_argument("--survey", choices=["ztf", "lsst"], default="lsst")
    ap.add_argument("--limit", type=int, default=100)
    ap.add_argument("--timeout", type=float, default=60.0)
    args = ap.parse_args()

    if args.check_creds:
        return 0 if check_creds() else 2
    if args.prune:
        prune(args.outdir, args.retention_days)
        return 0
    if args.poll_livestream:
        return poll_livestream(args.survey, args.limit, args.outdir, args.timeout)
    if args.topic:
        return pull(args.topic, args.outdir)
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
