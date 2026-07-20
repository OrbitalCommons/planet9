# 05 — Layer 4: slow-mover search (Phase 4)

> **Status 2026-07-19: implemented** (`p9_rubin_watch::linker` +
> `export_linker_input.py` + `link`/`--calibrate` CLI; RUNBOOK-linker.md).
> One algorithmic refinement vs the spec below, learned from testing: the
> 5-parameter refit searches the FULL 1/d range (coarse 60-point scan +
> fine refinement), because residual-parallax chains link a true object at
> badly wrong hypotheses too — the cluster's hypothesis is a proposal, not
> a prior. SR-6 gate passes as an offline test (≥ 90% for V ≤ 23, ≥ 3
> nights); live smoke on 4,509 selection-A detections (DE421) yields 0
> candidates as single-night physics demands. Awaiting selection-B
> multi-night data for the first real search.

The custom science layer: link unassociated DIA sources across nights into
candidates moving like a body at d ∈ [300, 1000] AU. This is the regime
where Rubin's own tracklet-based expectations break down and where the
corrected search hull says the surviving Planet Nine probability lives
(V ≈ 20.5–23 at 400–900 AU).

## Target kinematics

Apparent motion of a distant body is dominated by Earth's reflex:

| d (AU) | rate at opposition | motion in a 33-min visit pair | own orbital motion |
|---|---|---|---|
| 300 | 11.2″/day | 0.26″ | ~0.7″/day |
| 600 | 5.7″/day | 0.13″ | ~0.24″/day |
| 1000 | 3.4″/day | 0.08″ | ~0.11″/day |

Consequences that shape everything below:

1. **No intra-night tracklet.** 0.08–0.26″ over a visit pair is within the
   astrometric scatter — a P9-class object is a *stationary transient*
   within a night. Any pipeline requiring measurable intra-night motion is
   blind to it by construction; ours must link *nights*, not exposures.
2. **Rate windows are epoch-dependent.** The reflex rate sweeps from the
   opposition maximum through **zero at the stationary points** (twice per
   year) — near which the object is indistinguishable from a static
   transient for days–weeks, and *closer* objects (a 100-AU body near its
   own stationary point) sweep *through* our 2–14″/day window. Rate cuts
   alone therefore neither guarantee completeness nor purity: the linker
   must fit the full parallactic model, not a linear rate (SR-4).
3. **Brightness is not the hard part.** Single-visit r ≈ 24.3 covers the
   whole residual band with margin; SR-5's V ≤ 23.5 requirement is set by
   vetting robustness, not detectability.

## Algorithm: distance-hypothesis linking

Heliolinc-style transform, reusing the workspace's ephemeris machinery
(`p9_core::coords::observer::EarthProvider`, `apparent_position*`,
`candidate_pair` — built for the IRAS↔AKARI epoch pairing, which is the
same physics at a 20-year baseline):

```
for each tile (healpix nside 8, ~53 deg²) and 90-day window W:
  D ← unassociated records in tile±border(0.5°) over W        (design/04)
  for each hypothesis 1/d on the grid (below):
    for each detection: invert topocentric→heliocentric direction
      assuming distance d (Earth state from EarthProvider at each mjd)
    cluster in (λ, β) with tolerance θ_link, allowing a linear drift
      (dλ/dt, dβ/dt) ≤ the max heliocentric rate at d      [5-param model:
      λ0, β0, 1/d, dλ/dt, dβ/dt]
    clusters with ≥ 3 distinct nights → least-squares refit of the 5
      parameters over the raw (ra, dec, mjd); keep χ²/dof ≤ 3
  merge candidates found under adjacent hypotheses (same members ⇒ one
  candidate; record the 1/d interval that links it)
```

**Hypothesis grid.** A wrong 1/d displaces the transformed positions by
≈ B⊥·Δ(1/d) (B⊥ = Earth's projected baseline over the window, ≤ 2 AU for
90 days). Requiring < 2″ (≈ 9.7 µrad) residual:

- 90-day window: Δ(1/d) ≤ 4.8×10⁻⁶ AU⁻¹ over the range 1/1000–1/300 =
  2.33×10⁻³ AU⁻¹ ⇒ **≈ 480 hypotheses**.
- 14-day window (B⊥ ≲ 0.5 AU): ⇒ ≈ 120 hypotheses.

Per (tile, hypothesis), clustering ~10³–10⁴ points is a sort + neighbor
scan — the full footprint per month is minutes of CPU, embarrassingly
parallel over tiles. No performance heroics needed; correctness and
bookkeeping dominate.

**Window/cadence choice.** Primary: rolling 90-day windows stepped monthly,
centered away from each tile's stationary epochs (computed from the tile's
ecliptic longitude); secondary short windows (14 d) near opposition for the
highest-rate close-in tail. Windows spanning a stationary point get a
purity flag: static contaminants fit d→∞ there.

## Candidate scoring and auto-vetting chain

Score = f(nights, χ², arc length, magnitude consistency). A candidate must
survive, in order:

1. **Static veto**: the 5-param fit must beat (Δχ² ≥ 25) a static-source
   fit; kills variables/dipoles near stationary epochs.
2. **Known-object veto**: ssObjectId (already absent by selection), then
   astcheck against current MPCORB, then SkyBoT cone at each epoch (all
   recorded on the candidate record).
3. **Photometric consistency**: band-corrected magnitudes consistent with
   one reflected-light source at the fitted d over the arc
   (`p9_core::analysis::photometry`; the fitted d and the mass–albedo band
   must overlap the physical envelope — a labeled plausibility cut, wide).
4. **Cutout inspection set**: fetch stamps from the broker for every member
   (the only pixel-level step), auto-reject obvious artifacts (edge,
   saturation flags), package the rest for human review.
5. **Archival precovery**: predicted positions at PS1/DES/ZTF/DECaLS epochs
   using the fitted elements; a bright-enough predicted precovery that is
   *absent* in a clean archival image is strong evidence against; a
   *present* one is a discovery-grade confirmation. Depths come from the
   shared `SURVEY_DEPTHS` table.

Human review budget: ≤ ~10 candidates/month (SR from design/00); if the
chain leaks more, tighten score thresholds before loosening anything else.

**No MPC submission** in this phase (non-goal; OQ-1 in design/09). Output
is `rubin_watch/candidates/<id>.json` + review issue.

## Injection / recovery calibration (SR-6)

Truth in advertising for any exclusion claim we ever make from this search:

- Draw synthetic P9s from the ensemble posterior (`p9-search-hull`'s
  `Draws` machinery — same clouds as the hull).
- Generate synthetic DIA records using the *real* per-tile visit history
  from the Layer-2 coverage map (mjd, band, depth) + astrometric noise +
  the per-visit logistic detection efficiency.
- Inject into real tiles (so contamination and window structure are real),
  run the full linker + auto-vetting, measure recovery vs (V, d, nights,
  elongation).
- Gate: ≥ 90% recovery for V ≤ 23 with ≥ 3 covered nights in-window;
  publish the full efficiency surface with the monthly report.

## Hazards (each carries a mitigation, none is waved off)

| Hazard | Effect | Mitigation |
|---|---|---|
| Template self-subtraction (year 2+): the mover enters the coadd template, each visit then shows a positive/negative **dipole** | flux underestimated; naive linking sees two sources | year-1 data is clean by construction; from template-era onward, extend intake to keep negative-flux DIA sources and teach the static veto the dipole signature (positive+negative pair offset consistent with fitted motion × template epoch span); coverage map's `template_epoch_flag` (design/03) marks affected pixels |
| Stationary points | P9 static for weeks (missed links); static transients fit as d→∞ movers (false links) | window placement by tile ecliptic longitude; purity flag; static-veto Δχ² |
| Closer objects sweeping the rate window near their stationary points | false candidates at wrong d | full 5-param fit separates them (their parallax curvature fits small d, outside grid) — verify with injected 80–150 AU decoys in the calibration suite |
| Astrometric systematics per detector/band | fake curvature | fit includes per-band zero-point nuisance only if residuals demand it; start without |
| Broker reprocessing / alert revisions | duplicate or shifted records | store dedup on alert_id; linker consumes latest per dia_source_id |
| SSP already found it | wasted review | the known-object veto includes *Rubin's own* MPC submissions (weekly MPCORB refresh) — by construction we only surface what SSP hasn't designated |

## Honest scope statement

Rubin SSP + the community will link essentially everything at d ≲ 100–200
AU. This layer's *unique* contribution is the 300–1000 AU, no-tracklet,
near-stationary regime over our priority footprint — plus a measured
completeness there, which converts "we didn't find it" into a defensible
exclusion layer for the hull. That statement, with the efficiency surface,
is the deliverable even if no candidate ever survives vetting.
