#!/usr/bin/env python3
"""Submit a Fink LSST Data Transfer job programmatically.

The portal (lsst.fink-portal.org/download) is a Plotly Dash app whose
submit flow is a single unauthenticated callback POST to
`/_dash-update-component`; this script drives exactly that callback with
the same state the browser would send (verified against the portal source,
astrolabsoftware/lsst.fink-portal.org apps/datatransfer.py). The response
carries the server-generated topic name (`ftransfer_lsst_...`), which is
then consumable with our registered Kafka credentials via
`fink_pull.py --topic`.

Block negation follows the portal's own encoding: the literal label
"NOT <block>" (two clicks in the UI). Content options: "Full packet",
"Medium packet", "Light static packet", "Light SSO packet", the LSST
nested sections, or individual fields.

Selection-B example (the unassociated residue for the slow-mover search):
  fink_submit_job.py --start 2026-07-16 --stop 2026-07-18 \
      --block "NOT b_is_solar_system" --content "Light static packet" \
      --extra-cond "diaSource.reliability > 0.9;"

Etiquette: one job per invocation, exactly what the UI would do; jobs run
on Fink's Spark cluster, so keep requests scoped (design/08).
"""

import argparse
import json
import sys
import urllib.request

PORTAL = "https://lsst.fink-portal.org"
STATS_API = "https://api.lsst.fink-portal.org/api/v1/statistics"


def archived_nights(start: str, stop: str) -> list[str]:
    """Nights in [start, stop] present in Fink's statistics table — the same
    gate the server-side Spark job applies per night (check_path_exist); a
    night absent here is silently dropped, and a range with NO archived night
    makes the job sys.exit(1) before creating any topic (the batch-8/20
    failure: the archive stalled at 2026-07-14 while alerts kept streaming)."""
    with urllib.request.urlopen(
        urllib.request.Request(
            STATS_API,
            data=json.dumps({"date": "", "columns": "f:alerts",
                             "output-format": "json"}).encode(),
            headers={"Content-Type": "application/json"},
        ),
        timeout=60,
    ) as r:
        nights = {row["f:night"] for row in json.load(r)}
    lo, hi = start.replace("-", ""), stop.replace("-", "")
    return sorted(n for n in nights if lo <= n <= hi)


def fetch_output_spec() -> tuple[str, list]:
    """Read /_dash-dependencies and return the exact multi-output spec of the
    submit callback (the notification property carries a build hash that must
    match verbatim)."""
    with urllib.request.urlopen(f"{PORTAL}/_dash-dependencies", timeout=30) as r:
        deps = json.load(r)
    for d in deps:
        ins = [f"{i['id']}.{i['property']}" for i in d.get("inputs", [])]
        if ins == ["submit_datatransfer.n_clicks"]:
            out = d["output"]
            # "..a.b...c.d.." -> [(id, prop), ...]
            parts = [p for p in out.strip(".").split("...") if p]
            outputs = [
                {"id": p.split(".")[0], "property": p.split(".", 1)[1]} for p in parts
            ]
            return out, outputs
    raise SystemExit("submit_datatransfer callback not found in dash dependencies")


def submit(start: str, stop: str, blocks: list, tags: list, content: list,
           extra_cond: str | None) -> dict:
    out_spec, outputs = fetch_output_spec()
    state = [
        {"id": "date-range-picker", "property": "value", "value": [start, stop]},
        {"id": "tag_select", "property": "data", "value": tags},
        {"id": "blocks_select", "property": "data", "value": blocks},
        {"id": "field_select", "property": "value", "value": content},
        {"id": "extra_cond", "property": "value", "value": extra_cond},
        {"id": "object-catalog", "property": "data", "value": None},
        {"id": "upload-data", "property": "filename", "value": None},
        {"id": "ra-column", "property": "value", "value": None},
        {"id": "dec-column", "property": "value", "value": None},
        {"id": "radius_xmatch", "property": "value", "value": None},
        {"id": "id-column", "property": "value", "value": None},
    ]
    payload = {
        "output": out_spec,
        "outputs": outputs,
        "inputs": [{"id": "submit_datatransfer", "property": "n_clicks", "value": 1}],
        "changedPropIds": ["submit_datatransfer.n_clicks"],
        "state": state,
    }
    req = urllib.request.Request(
        f"{PORTAL}/_dash-update-component",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", required=True, help="YYYY-MM-DD (both dates required; "
                    "a missing stop date crashes the server-side date parser)")
    ap.add_argument("--stop", required=True, help="YYYY-MM-DD")
    ap.add_argument("--block", action="append", default=[],
                    help='e.g. "b_is_solar_system" or "NOT b_is_solar_system"; repeatable')
    ap.add_argument("--tag", action="append", default=[],
                    help="Fink filter tag; repeatable")
    ap.add_argument("--content", action="append", default=[],
                    help='e.g. "Light SSO packet", "Light static packet"; '
                    "default server-side: Full packet")
    ap.add_argument("--extra-cond", default=None,
                    help='semicolon-terminated SQL, e.g. "diaSource.reliability > 0.9;"')
    args = ap.parse_args()

    nights = archived_nights(args.start, args.stop)
    if not nights:
        print(f"no archived nights in [{args.start}, {args.stop}] — the job "
              "would exit(1) server-side without creating a topic; pick nights "
              f"listed by {STATS_API}", file=sys.stderr)
        return 2
    print(f"archived nights in range: {' '.join(nights)}", file=sys.stderr)

    resp = submit(args.start, args.stop, args.block, args.tag, args.content,
                  args.extra_cond)
    r = resp.get("response", {})
    topic = r.get("topic_name", {}).get("children")
    batch = r.get("batch_id", {}).get("children")
    notes = r.get("notification-container", {})
    print(json.dumps({"topic": topic, "batch_id": batch}, indent=1))
    if notes:
        for n in next(iter(notes.values()), []) if isinstance(notes, dict) else []:
            msg = n.get("message") if isinstance(n, dict) else n
            print(f"note: {msg if isinstance(msg, str) else '(rich message)'}",
                  file=sys.stderr)
    if not topic:
        print("no topic returned — inspect full response:", file=sys.stderr)
        print(json.dumps(resp, indent=1)[:2000], file=sys.stderr)
        return 2
    print(f"\nconsume with:\n  fink_pull.py --topic {topic}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
