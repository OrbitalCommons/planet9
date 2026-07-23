# Runbook — Fink broker access (Layer 3)

Registration is done (2026-07-17, fink-client 11.0). Credentials are
low-sensitivity (PLAINTEXT Kafka, no password) but keep the email in
1Password; the working copies live in `~/.finkclient/{lsst,ztf}_credentials.yml`
and are NOT in the repo.

## Environment

```bash
uv venv ~/.venvs/fink && uv pip install --python ~/.venvs/fink/bin/python "fink-client>=11"
```

(System python is PEP-668 managed and `python3 -m venv` lacks ensurepip on
this box — use uv.)

## Re-registering / changing topics

```bash
~/.venvs/fink/bin/fink_client_register -survey lsst \
  -username matthew -group_id matthew_cosmicfrontier \
  -servers kafka-lsst.fink-broker.org:24499 \
  -mytopics fink_uniform_sample_lsst fink_hostless_candidate_lsst -maxtimeout 30
# same for -survey ztf with kafka-ztf.fink-broker.org:24499
```

Topic choices can change any time (re-register). Discover what exists with
Kafka metadata (no docs needed):

```bash
~/.venvs/fink/bin/python -c "
from confluent_kafka import Consumer
c = Consumer({'bootstrap.servers': 'kafka-lsst.fink-broker.org:24499',
              'group.id': 'matthew_cosmicfrontier'})
print(sorted(t for t in c.list_topics(timeout=15).topics
             if t.startswith('fink_')))"
```

## Topic map (as verified 2026-07-17)

- **LSST curated (11 topics)**: extragalactic/SN oriented plus
  `fink_uniform_sample_lsst` and `fink_hostless_candidate_lsst`.
  **No LSST SSO livestream topic exists yet** — the LSST solar-system path
  is Data Transfer (`b_is_solar_system` block / Light-SSO packets), per
  design/04. `fink_hostless_candidate_lsst` (no host crossmatch) partially
  overlaps our selection-B and is worth watching: a P9-class stationary
  transient would qualify.
- **ZTF**: `fink_sso_ztf_candidates_ztf` (SSO-shaped movers) and
  `fink_sso_fink_candidates_ztf` (Fink's own candidates, i.e. unassociated)
  are live and flowing — the end-to-end testbed for our consumer code.

## Polling (livestream)

```bash
timeout 120 ~/.venvs/fink/bin/fink_consumer -survey lsst --display -limit 5
timeout 120 ~/.venvs/fink/bin/fink_consumer -survey ztf  --display -limit 5
# --save -outdir DIR to keep Avro packets
```

Queue retention is ~7 days; streams are replayable within it. Verified
2026-07-17: live LSST hostless alerts (diaObjectId-keyed) and ZTF SSO
candidates both delivered.

## Landing livestream alerts into the store

```bash
~/.venvs/fink/bin/python scripts/rubin_watch/fink_pull.py   --poll-livestream --survey ztf --limit 100 --timeout 60
```

Writes `night=YYYYMMDD/tile=<nside8>/part-NNNN.parquet` per design/07
(dedup on read via alert_id). Verified 2026-07-17: 12 live ZTF SSO alerts
landed across two nights, healpix-tiled, full column set, known-object
associations (`ss_name`) intact.

## Data Transfer (the LSST SSO path)

Job creation is fully scripted — the portal's submit flow is an
unauthenticated Dash callback, driven directly by
`scripts/rubin_watch/fink_submit_job.py` (verified against the portal
source; first programmatic job: batch 8, 2026-07-20):

```bash
# selection A (SSO-tagged)
~/.venvs/fink/bin/python scripts/rubin_watch/fink_submit_job.py   --start 2026-07-16 --stop 2026-07-18   --block b_is_solar_system --content "Light SSO packet"

# selection B (unassociated residue, the slow-mover feed)
~/.venvs/fink/bin/python scripts/rubin_watch/fink_submit_job.py   --start 2026-07-16 --stop 2026-07-18   --block "NOT b_is_solar_system" --content "Light static packet"   --extra-cond "diaSource.reliability > 0.9;"
```

Both dates are REQUIRED (a missing stop date crashes the server-side date
parser — the Batch-5 failure). Negation is the literal label
`NOT <block>`. One job per invocation; their Spark cluster does the work,
so keep requests scoped.

**Nights must exist in Fink's archive.** The Spark job gates every night on
the statistics API (`api.lsst.fink-portal.org/api/v1/statistics`); a night
missing there is silently dropped, and if the whole range is missing the
job `sys.exit(1)`s before creating any topic — no error surfaces anywhere
(batches 8/20: the archive stalled at 2026-07-14 while livestream alerts
kept flowing, so "recent" nights were unfetchable). The submitter now
pre-checks the range and prints the archived nights; if it reports none,
pick different dates — resubmitting the same range only wastes their
cluster. Archive lag vs livestream is a known divergence; if it persists
more than a week, mail contact@fink-broker.org.

1. (Manual fallback: https://lsst.fink-portal.org/download with the same
   choices.)
2. The job yields a private topic `ftransfer_lsst_<date>_<id>` (lives ~7
   days; allow minutes for the Spark job to start streaming). Consume with:

```bash
~/.venvs/fink/bin/fink_datatransfer -survey lsst \
  -topic ftransfer_lsst_... -outdir ~/.cache/p9-rubin-watch/alerts/raw \
  -nconsumers 4 --dump_schemas --verbose
```

3. Re-partition to the design/07 parquet layout with
   `scripts/rubin_watch/fink_pull.py` (volume alarm enforced there).

Portal job creation is a browser step (login required); everything after
the topic name is scriptable.
