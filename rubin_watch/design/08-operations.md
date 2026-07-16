# 08 — Operations

## Cadence

| Task | Cadence | Runtime | Command (target CLI) |
|---|---|---|---|
| MPC sync + battery | weekly (Mon) | ~1 min + battery | `cargo run -p p9-rubin-watch -- mpc-sync` |
| Frozen-table drift (piggybacked) | with mpc-sync | ~30 s | (included) |
| Fink batch pull | weekly → nightly once stable | ~min, network-bound | `python scripts/rubin_watch/fink_pull.py --nights …` |
| Coverage rebuild + hull/viability regen | monthly (1st) | ~5 min | `cargo run -p p9-rubin-watch -- coverage && cargo run --release -p p9-search-hull && cargo run --release -p p9-viability` |
| Slow-mover link + vetting | monthly, after coverage | ~10 min CPU | `cargo run --release -p p9-rubin-watch -- link --window 90d` |
| Injection/recovery calibration | monthly, with link | ~30 min | `… -- link --calibrate` |

Manual-first: every task is a runbook command producing a branch + PR, in
the paper-watch style. Cron only after ≥ 3 clean manual cycles per task,
and cron still opens PRs — nothing lands on main without the gate.

## Git and PR conventions

- Branch per run: `meawoppl/rubin-watch-YYYY-MM-DD[-layer]`.
- One PR per run; body carries the delta-report summary; `Fixes #N` only
  when a run closes a tracked issue.
- Ledger merge policy: the ledger is regenerated, not hand-merged — on
  conflict, re-run the sync on a fresh branch (append-only history +
  canonical serialization make this safe and cheap).
- Gate: the standard full-suite worktree gate applies to every PR,
  including docs/ledger-only runs (fixture tests exercise the parsers, so
  a format drift is caught by the gate, not by the next live run).
- Reports and committed artifacts follow the repo attribution rule: no
  AI attribution anywhere.

## Credentials

- Fink username/group (PLAINTEXT Kafka — treat as low-sensitivity but do
  not commit): stored in 1Password, injected via `op inject` into env
  (`FINK_USERNAME`, `FINK_GROUP`) per the global CLAUDE.md pattern. If the
  1Password app connection fails, the runbook says: open the 1Password app
  first.
- No other layer holds credentials (MPC, SBDB, SkyBoT are anonymous).

## Storage budget

| Store | Location | Budget | Enforcement |
|---|---|---|---|
| Ledger + coverage + candidates + reports | repo (committed) | ≤ ~5 MB total; coverage arrays dominate (~1 MB) | review at PR |
| Alert parquet | `~/.cache/p9-rubin-watch/alerts/` | ≤ 10 GB (120-night hot window × ≤ 50 MB) | `fink_pull.py --prune` drops partitions beyond retention |
| Cutouts | `~/.cache/p9-rubin-watch/cutouts/` | ≤ 1 GB | fetched per-candidate only; pruned with candidate resolution |

## Failure and alerting policy

- Any task exiting `2` (source failure, design/02) retries next scheduled
  run; **3 consecutive** failures of the same task ⇒ open an ops issue
  labeled `rubin-watch-ops`.
- Volume alarm: a Fink pull > 3× the nightly budget aborts before write
  and opens an ops issue (filter regression, not data).
- Schema alarm: any `schema_version` mismatch or parquet column drift
  aborts loudly; never coerce.
- Degraded runs (`degraded:` nonempty) are allowed to land but the report
  header says so; two consecutive degraded runs from the same source ⇒
  ops issue.
- The linker never blocks the MPC layer: layers fail independently.

## External-service etiquette

- MPC file: one conditional GET (If-Modified-Since) per run.
- SBDB: ≤ 1 request/s (existing SbdbClient behavior), single-object
  lookups only for new/refresh-due objects.
- SkyBoT: batch cone queries, ≤ 1/s, only for linker candidates (tens per
  month) — never bulk sweeps (astcheck against the local MPCORB handles
  bulk).
- Fink: use Data Transfer batch jobs as designed; livestream group-id
  registered once; respect their published limits when they publish them.
- Attribution: reports and any publication acknowledge Rubin/LSST alerts,
  Fink, MPC, IMCCE SkyBoT per each service's citation policy (tracked in
  the report template).

## Runbook index (to be written with each phase's implementation)

- `RUNBOOK-mpc.md` — weekly sync: command, expected output, how to read
  the delta report, frozen-table adoption procedure (human PR).
- `RUNBOOK-coverage.md` — monthly rebuild + hull/viability regen +
  headline-number update procedure.
- `RUNBOOK-fink.md` — registration, credential injection, pull, prune,
  volume checks.
- `RUNBOOK-linker.md` — monthly run, calibration read-out, candidate
  review checklist (what a human checks before `human_review` →
  `confirmed`), escalation path (OQ-1 before any external report).
