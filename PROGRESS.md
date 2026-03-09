# Planet 9 — Implementation Progress

## Overview

Re-implementation of all numerical models from Batygin, Brown et al. Planet Nine papers in Rust, using the `starfield` crate for astronomical primitives and a custom N-body integration core.

## Workspace Structure

| Crate | Paper | Status |
|-------|-------|--------|
| `p9-core` | Shared integrator, orbital mechanics, I/O | MERGED |
| `p9-2016-evidence` | Evidence for a Distant Giant Planet (Batygin & Brown 2016) | MERGED |
| `p9-2016-constraints` | Observational Constraints on Planet Nine (Brown & Batygin 2016) | MERGED |
| `p9-2016-obliquity` | Solar Obliquity Induced by Planet Nine (Bailey+ 2016) | IN PROGRESS |
| `p9-2016-inclined-tnos` | Highly Inclined TNOs by Planet Nine (Batygin & Brown 2016) | NOT STARTED |
| `p9-2017-bias` | Observational Bias and Clustering (Brown 2017) | NOT STARTED |
| `p9-2018-kuiper-belt` | Generation of Distant Kuiper Belt (Khain+ 2018) | NOT STARTED |
| `p9-2018-resonance` | Resonance-based Planet Nine Search (Bailey+ 2018) | NOT STARTED |
| `p9-2019-clustering` | Orbital Clustering in the Distant Solar System (Brown & Batygin 2019) | NOT STARTED |
| `p9-2019-review` | The Planet Nine Hypothesis Review (Batygin+ 2019) | NOT STARTED |
| `p9-2021-oort-cloud` | Injection of Inner Oort Cloud (Batygin & Brown 2021) | NOT STARTED |
| `p9-2021-orbit` | The Orbit of Planet Nine (Brown & Batygin 2021) | NOT STARTED |
| `p9-2021-ztf` | ZTF Search for Planet Nine (Brown & Batygin 2021) | NOT STARTED |
| `p9-2022-des` | DES Limits on Planet Nine (Belyakov+ 2022) | NOT STARTED |
| `p9-2024-panstarrs` | Pan-STARRS1 Search (Brown+ 2024) | NOT STARTED |
| `p9-2024-neptune-crossing` | Neptune-Crossing TNOs (Batygin+ 2024) | NOT STARTED |

## Progress Log

### Iteration 0 — Project Setup
- [x] Downloaded 15 Planet Nine papers (PDF + TeX source) from arXiv
- [x] Created detailed summaries of all papers (SUMMARIES_2016.md, SUMMARIES_2017_2019.md, SUMMARIES_2021_2024.md)
- [x] Identified all numerical models requiring re-implementation
- [x] Designed workspace structure and implementation plan
- [x] Initialize Cargo workspace with all subcrates
- [x] Set up `p9-core` with integrator scaffolding

### Iteration 1 — p9-core Implementation (2026-03-09)
- [x] Constants: GM values (DE440), J2/J4, unit conversions, galactic tidal density
- [x] Types: StateVector, OrbitalElements, MassiveBody, Particle, P9Params, SimConfig
- [x] Coordinate conversions: elements↔Cartesian (Murray & Dermott), Kepler equation solver
- [x] Democratic heliocentric coordinates (WHM requirement)
- [x] Kepler drift: universal variable formulation with Stumpff functions
- [x] Perturbation kicks: direct+indirect, J2/J4 oblateness, parallel (rayon)
- [x] Wisdom-Holman symplectic integrator (kick-drift-kick leapfrog)
- [x] Bulirsch-Stoer adaptive integrator (modified midpoint + Richardson extrapolation)
- [x] Hybrid integrator (WHM/BS switching on Hill sphere proximity)
- [x] J2 secular: orbit-averaged giant planet quadrupole
- [x] Galactic tide: vertical Milky Way disk acceleration
- [x] Initial conditions: giant planet state vectors (J2000), scattered disk generators
- [x] Analysis: orbital element snapshots, longitude of perihelion, Δϖ computation
- [x] Secular Hamiltonian: quadrupole (coplanar + 3D), phase-space portrait generation
- [x] Visualization: SVG heatmap + scatter plot generation, viridis colormap
- [x] All 14 unit tests passing
- [x] PR #17 merged — 2026-03-09

### Iteration 2 — p9-2016-evidence: Batygin & Brown (2016) (2026-03-09)
- [x] KBO elements: hardcoded orbital elements for 6 stable KBOs, clone generation
- [x] Octupole Hamiltonian: coplanar + 3D extensions to p9-core's quadrupole
- [x] Phase-space portraits: N-body trajectory tracing, alignment classification
- [x] Scattered disk simulation: planar + 3D configs, snapshot recording, clustering statistics
- [x] Resonance detection: MMR identification (11 known resonances), libration analysis
- [x] SVG plots: Figure 1 (KBO clustering), phase portraits, scattered disk clustering
- [x] All 30 unit tests passing (16 evidence + 14 core)
- [x] PR #18 merged — 2026-03-09

### Iteration 3 — p9-2016-constraints: Brown & Batygin (2016) (2026-03-09)
- [x] Parameter grid: 19×9×5 (a, e, mass) survey with perihelion filtering
- [x] Acceptance criteria: n_survivors ≥ 7, clustering > 0.5, high-q > 0.05
- [x] Clustering metrics: Rayleigh test, anti-aligned fraction, confinement probability
- [x] Inclination survey: 8×12 (i₉, ω₉) grid, pole direction, angular separation
- [x] Detection limits: V magnitude from mass/albedo/distance, mass-radius relation
- [x] Sky position: ecliptic coordinate computation along orbit
- [x] Brightness curves: V mag vs true anomaly with survey depth overlays
- [x] SVG plots: acceptance heatmap (Fig 2), brightness curve (Fig 5)
- [x] All 47 unit tests passing (17 constraints + 16 evidence + 14 core)
- [x] PR #19 merged — 2026-03-09

### Iteration 4 — p9-2016-obliquity: Bailey, Batygin & Brown (2016) (2026-03-09)
- [x] Solar model: n=3 polytrope, I_hat=0.08, k₂=0.01, Skumanich spin-down
- [x] Vector secular Hamiltonian: 3-body quadrupole (gp + P9 + solar spin)
- [x] Angular momentum coupling: time-dependent solar ring, effective masses
- [x] RK4 integration with spin-down-driven coupling evolution
- [x] Parameter survey: bisection search for required i₉ to produce 6° obliquity
- [x] SVG plots: required inclination contours (Fig 2), obliquity evolution
- [x] All 60 unit tests passing (13 obliquity + 17 constraints + 16 evidence + 14 core)
