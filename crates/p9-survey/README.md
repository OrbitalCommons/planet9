# p9-survey

A "survey paper" in code: it gathers the Planet Nine orbit solutions reproduced
across this workspace and asks two questions numerically —

1. **Where is it most likely now?** For each solution it Monte-Carlos the
   unknown current orbital phase (drawn uniformly in mean anomaly → dwell-time
   weighted, so the planet is most often near aphelion) plus the solution's
   assumed 1σ element spreads, propagates with `p9_core::elements_to_cartesian`,
   and bins the resulting RA/Dec into a sky probability map. Brightness is the
   reflected-light V from `p9_core::analysis::photometry`.
2. **Can a survey catch it?** Each telescope (`telescope.rs`) has a footprint
   (declination band + galactic-plane cut + coverage fraction) and a limiting
   magnitude; detection probability = P(in footprint **and** brighter than the
   depth), integrated over the same Monte Carlo.

Surveyed solutions (all sourced from reproduction crates, not re-typed):

| Solution | source |
|---|---|
| 2016 Batygin & Brown (nominal + inclined-TNO variant) | `p9_core::types::P9Params` |
| 2019 Batygin et al. (review best-fit) | `p9_core::types::P9Params` |
| 2021 Brown & Batygin (MCMC median) | `p9_core::types::P9Params` |
| 2024 Siraj, Chyba & Tremaine (independent) | `p9_2024_siraj_orbit::reference` |

## Regenerate the figure

```sh
# 1. Rust computes everything and writes the typed dataset:
cargo run -p p9-survey --bin survey -- \
    --out crates/p9-survey/p9_survey_data.json --samples 300000 --seed 2026

# 2. Python renders it (no astronomy of its own — only what Rust emitted):
python3 scripts/plot_survey.py \
    crates/p9-survey/p9_survey_data.json crates/p9-survey/p9_survey_sky.svg
```

The JSON is a `schema::SurveyDataset` (see `src/schema.rs`); the plotter mirrors
that exact contract with typed dataclasses and rejects any file whose
`schema_version`, keys, or field types don't match. Bump `SCHEMA_VERSION` on
**both** sides for any breaking change. The committed `p9_survey_sky.svg` is the
rendered output; the large `p9_survey_data.json` is reproducible and not checked
in.

The brightness here is reflected-light optical/NIR; thermal-IR detectability of
a *cold* Planet Nine is treated (and shown negligible at WISE W1) in
`p9-2018-wise-search`.
