# p9-search-hull

Where to point next. This crate crosses the Planet Nine **current-position
posterior** against the combined **coverage hull** of every wide optical search,
and reports the *residual*: the high-probability sky that has **not** yet been
searched deeply enough — restricted to what a small all-sky space telescope
(the JBT 0.5 m / SPENCER concept, V ≈ 24.5) could actually detect.

It reuses the workspace's own models — `p9-core` photometry and coordinates,
and the orbit-solution catalog + sky grid from `p9-survey` — so the map tracks
the reproductions, not a hand-drawn cartoon.

## The model

- **Sizes → magnitudes.** Mass → Neptune-anchored radius → absolute `H` →
  apparent reflected `V` (`p9_core::analysis::photometry`).
- **Directions + distances.** Orbital elements → heliocentric position →
  equatorial RA/Dec + heliocentric distance (`p9_core::coords`).
- **Where it could be.** An equal-weight **ensemble** over every orbit solution
  in `p9_survey::studies`, each Monte-Carlo sampled (1σ element jitter, uniform
  mean anomaly = dwell-weighted "where is it now"). Per sky cell we keep the
  probability and the probability-weighted predicted distance and `V`.
- **What's been observed.** Each wide optical survey (CRTS, ZTF, PS1 3π, DES) as
  a footprint (declination band + galactic cut + areal coverage) and a depth
  from `p9_core::analysis::surveys::SURVEY_DEPTHS`. Combined per direction into
  the deepest search that reaches it — the coverage hull — and a per-survey-OR
  exclusion probability for a given predicted `V`.
- **Residual + ranking.** `residual = posterior × (1 − P_excluded)`, then ranked
  over cells our space telescope can both see (footprint) and detect (`V ≤
  24.5`).
- **Parameter-space hull.** The (distance × apparent `V`) plane with the
  posterior cloud and the survey-depth cuts that bound what has been seen.

## Run

```
cargo run -p p9-search-hull            # writes figures/search_hull.json
PYTHONNOUSERSITE=1 python3 scripts/plot_search_hull.py   # writes the two SVGs
```

## Headline (current inputs)

- **≈65%** of the Planet Nine posterior is already excluded by past optical
  surveys; **≈35%** remains un-searched.
- **≈29%** of the total posterior is both un-searched **and** reachable at
  V ≤ 24.5 from space.
- The deepest all-sky survey (PS1, r ≈ 21.5) reaches ~600 AU for a 6.2 M⊕
  planet; the V ≈ 24.5 space telescope reaches ~1200 AU.
- Top pointings cluster near **RA ≈ 55–75°, Dec ≈ +8–14°** (V ≈ 21.6,
  d ≈ 670 AU) — the predicted-aphelion arc, just past PS1's depth.

The JSON also carries the **per-study** prediction clouds (each orbit solution's
thinned position + distance + V draws) and a **Rubin/LSST** reach layer
(footprint + depth + per-cell reachability), so the figures can show each study
on its own and compare the ground baseline.

Outputs: `figures/search_hull.json`, `figures/p9_search_sky.svg` (coverage hull +
where-to-look, with the LSST footprint), `figures/p9_reach_hull.svg`
(distance × V, per-study clouds), `figures/p9_study_clouds.svg` (per-study
clouds on the sky and in distance).
