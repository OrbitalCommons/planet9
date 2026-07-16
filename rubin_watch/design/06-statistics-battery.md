# 06 — Statistics battery

What actually recomputes when Layer 1 reports a delta, how results are
reported, and which downstream artifacts regenerate. The battery reuses the
paper crates as libraries — no statistic is reimplemented in the watcher.

## Battery table

Run order is fixed; each row states its input sample rule and the crate
function(s) invoked. Every clustering statistic is computed twice — vetted
tier only, and vetted+provisional — and reported side by side (SR-2).

| # | Statistic | Sample rule | Machinery |
|---|---|---|---|
| B1 | ϖ mean-resultant R̄, von Mises κ, Gaussian-equivalent σ | `etno_vetted_class` | `p9-2026-apsidal-clustering` estimator pipeline (`kappa_from_r_bar`, `p_value_to_sigma` via core stats) |
| B2 | Rayleigh + Kuiper p on ϖ, Ω, ω | same | `p9_core::analysis::circular` |
| B3 | Selection-conditioned clustering p under all three labeled stand-ins | same | `p9_core::analysis::selection::VarpiSelection` MC (the verdict-sensitivity machinery); report the spread, never a single number (SR-3) |
| B4 | Perihelion-gap continuous-null p, paper-epoch and current-epoch samples | `gap_band` + `gap_context` per p9-2021-perihelion-gap rules | `continuous_null_dip_p_value` |
| B5 | Neptune-crossing sample delta | `neptune_crossing_class` | membership check against `p9-2024-neptune-crossing::observed_tnos`; if changed, re-run ζ/KS headline |
| B6 | High-i / retrograde census delta | `high_inclination`, `retrograde` | table diff vs `p9-2016-inclined-tnos::known_objects` |
| B7 | Frozen-table drift | all | existing sbdb-refresh guard, all five guarded tables |
| B8 | Orbit-posterior staleness check | `etno_vetted_class` | rule-based flag, not a re-fit (below) |

**B8 rule.** The published orbit solutions (`p9_survey::studies`) were fit
to specific samples; we do not re-run anyone's MCMC automatically. Flag
`posterior_stale` when the vetted sample differs from `BROWN_2017_SAMPLE`
by ≥ 3 objects (added or removed) — the delta report then recommends
opening a re-derivation issue, and the hull's ensemble caveat is quoted in
any hull regeneration until resolved.

## Delta report

`reports/rubin-watch/YYYY-MM-DD-<layer>.md`, frontmatter per design/07.
Body structure (fixed, diffable):

1. **Trigger** — new/promoted/reclassified objects with elements, tier,
   tags, provenance.
2. **Battery results** — one table per B-row: previous value → new value
   (vetted) and (vetted+provisional), with the with/without-new-object
   column when the trigger is a single object.
3. **Interpretation guardrails** — auto-inserted boilerplate: selection
   spread (B3) restated; provisional-tier caveat; a–e degeneracy note for
   any 1-opposition object.
4. **Actions** — artifacts regenerated (below), issues opened, frozen-table
   adoption recommendations (human PR required).

## Regeneration decision rules

| Condition | Action |
|---|---|
| Any B1–B4 change ≥ 0.1σ equivalent (vetted tier) | regenerate `figures/search_hull.json` + `figures/viability.json`; note headline movement in report |
| `posterior_stale` flag newly set | issue: "re-derive orbit-solution inputs"; hull regen proceeds with caveat quoted |
| Coverage-map monthly rebuild moves excluded% ≥ 0.5 pt | hull + viability regen; report |
| Any frozen-table drift (B7) | drift section in report; fix-forward PR per the sbdb-refresh conventions |
| review-article scope counters change (paper crate count/year span) | nothing — `p9-review-article` derives them at build time |

## Worked example (to be executed as the Phase-1 inaugural report)

Trigger: 2025 LS2 (a = 523 AU, e = 0.917, q = 43.2 AU, i = 12.6°,
Ω = 144.0°, ω = 192.0°, ϖ = 336.0°; arc 350 d, soln 2 ⇒ tier
`provisional`, tag `etno_vetted_class`).

- B1/B2 expectation (qualitative, to be quantified by the run): ϖ = 336°
  sits ≈ 76° from the cluster mean ≈ 52° (≈ 2.2 circular-σ given σ ≈ 34°);
  including it in the provisional column lowers R̄ and κ and raises the
  Rayleigh p. The vetted column is unchanged (1 opposition).
- B3: the movement will differ by stand-in — under the cluster-aligned lobe
  an off-cluster object is *expected* (selection suppresses nothing at
  336°... report the actual three numbers); this is exactly the case the
  spread-reporting rule exists for.
- B4: q = 43.2 is below the 50–75 AU gap band → `gap_context` only; gap p
  unchanged.
- B8: |Δsample| = 1 → no staleness flag.

The point of the worked example in this document is the *shape* of the
report; the numbers come from the battery, never from prose.
