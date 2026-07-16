# Runbook — MPC watch (Layer 1)

Weekly sync of the distant-object census into `rubin_watch/etno_ledger.json`
with a statistics-battery delta report. Design: `design/02`, `design/06`.

## The run

```bash
git checkout -b meawoppl/rubin-watch-$(date +%F) main
cargo run --release -p p9-rubin-watch -- mpc-sync
```

Exit codes: `0` no changes (nothing to commit), `1` delta emitted (ledger +
report written), `2` source failure (nothing mutated; retry next week, ops
issue after 3 consecutive failures).

Dry-run (prints the report, touches nothing): add `--dry-run`.
Offline replay against fixtures: `--fixtures crates/p9-rubin-watch/tests/fixtures`.

## Reading the delta report (`reports/rubin-watch/YYYY-MM-DD-mpc.md`)

- **Trigger** — new/promoted/reclassified objects with elements and quality.
  A `discrepant` entry means MPC and SBDB disagree beyond the drift
  tolerances; the ledger keeps both element sets (`elements` = SBDB,
  `elements_alt` = MPC). Common for fresh Rubin discoveries riding the a–e
  degeneracy; it resolves on later oppositions.
- **B1/B2 table** — read the *vetted* row for anything you'd quote; the
  provisional column is advisory. The frozen-baseline row should be stable
  run over run (it is computed from the paper crate's own sample).
- **B3 spread** — if the three numbers straddle 0.05, the clustering verdict
  is currently selection-model-bound; say so wherever the number travels.
- **B4/B5/B6 lists** — frozen-table adoption candidates. Adoption is a
  human PR: update the relevant `p9_core::data::*` / paper-crate table with
  provenance comments, re-pin downstream tests, cite the ledger id.
- **B8 stale** — open (or link) a re-derivation issue for the orbit-solution
  inputs; hull regenerations quote the caveat until resolved.

## Landing it

```bash
git add rubin_watch/etno_ledger.json reports/rubin-watch/
git commit -m "Rubin watch MPC sync YYYY-MM-DD"
git push -u origin <branch> && gh pr create ...
```

Standard full-suite gate before merge. Ledger conflicts: never hand-merge —
re-run the sync on a fresh branch from main.
