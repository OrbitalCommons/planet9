#!/usr/bin/env python3
"""Fink Data Transfer pull for the Rubin watch (Layer 3, design/04).

Credential-gated: does nothing useful until Fink registration lands
(`rubin_watch/design/04-layer3-broker-intake.md`). Registration writes
~/.finkclient/lsst_credentials.yml via `fink_client_register -survey lsst ...`;
this script refuses to run without it rather than half-working.

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


def pull(topic: str, outdir: pathlib.Path) -> int:
    if not check_creds():
        return 2
    try:
        # Imported lazily: fink-client is only needed on the machine that pulls.
        from fink_client.consumer import AlertConsumer  # noqa: F401
    except ImportError:
        print("pip install fink-client (>=10, LSST-era) first", file=sys.stderr)
        return 2
    # Consumption loop lands with the Phase-3 spike once credentials exist:
    # fink_datatransfer is the reference consumer; this script will either
    # shell out to it or use AlertConsumer directly, then re-partition to the
    # design/07 parquet layout and enforce NIGHTLY_RECORD_BUDGET * ALARM_FACTOR.
    print(
        f"credentials ok; Phase-3 consumption loop not yet implemented "
        f"(topic {topic}, outdir {outdir}) — see design/04 acceptance tests",
        file=sys.stderr,
    )
    return 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-creds", action="store_true")
    ap.add_argument("--topic")
    ap.add_argument("--outdir", type=pathlib.Path,
                    default=pathlib.Path.home() / ".cache" / "p9-rubin-watch" / "alerts")
    ap.add_argument("--prune", action="store_true")
    ap.add_argument("--retention-days", type=int, default=120)
    args = ap.parse_args()

    if args.check_creds:
        return 0 if check_creds() else 2
    if args.prune:
        prune(args.outdir, args.retention_days)
        return 0
    if args.topic:
        return pull(args.topic, args.outdir)
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
