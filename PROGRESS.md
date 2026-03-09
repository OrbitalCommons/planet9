# Planet 9 — Implementation Progress

## Overview

Re-implementation of all numerical models from Batygin, Brown et al. Planet Nine papers in Rust, using the `starfield` crate for astronomical primitives and a custom N-body integration core.

## Workspace Structure

| Crate | Paper | Status |
|-------|-------|--------|
| `p9-core` | Shared integrator, orbital mechanics, I/O | IN PROGRESS |
| `p9-2016-evidence` | Evidence for a Distant Giant Planet (Batygin & Brown 2016) | NOT STARTED |
| `p9-2016-constraints` | Observational Constraints on Planet Nine (Brown & Batygin 2016) | NOT STARTED |
| `p9-2016-obliquity` | Solar Obliquity Induced by Planet Nine (Bailey+ 2016) | NOT STARTED |
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
