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

## Data Transfer (the LSST SSO path)

1. Create the job in the portal: https://lsst.fink-portal.org/download —
   pick nights, apply `b_is_solar_system` (selection A) or negate it +
   quality cuts (selection B), choose the Light-SSO packet schema.
2. The job yields a private topic `ftransfer_lsst_<date>_<id>` (lives ~7
   days). Consume with:

```bash
~/.venvs/fink/bin/fink_datatransfer -survey lsst \
  -topic ftransfer_lsst_... -outdir ~/.cache/p9-rubin-watch/alerts/raw \
  -nconsumers 4 --dump_schemas --verbose
```

3. Re-partition to the design/07 parquet layout with
   `scripts/rubin_watch/fink_pull.py` (volume alarm enforced there).

Portal job creation is a browser step (login required); everything after
the topic name is scriptable.
