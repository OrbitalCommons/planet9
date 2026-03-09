# Planet 9 — Rust Re-Implementation Plan

Re-implementation of every numerical model from the Batygin/Brown Planet Nine paper series. Each paper maps to a subcrate in a Cargo workspace. The `starfield` crate provides JPL ephemeris, coordinate transforms, and orbital element conversions. A shared `p9-core` crate provides the N-body integrator and common types.

---

## Shared Core: `p9-core`

### Components

**N-body Integrator** (mercury6 equivalent)
- Wisdom-Holman symplectic mapping with democratic heliocentric coordinates
- Bulirsch-Stoer adaptive integrator for close encounters
- Hybrid switching logic (Hill sphere proximity trigger)
- Configurable timestep (100 days to 1/20 Jupiter period)
- Accuracy parameter epsilon (10^-11 to 10^-12)

**Orbital Mechanics**
- Keplerian elements <-> Cartesian state vector conversion (leverage `starfield`)
- Secular perturbation Hamiltonian (quadrupole + octupole order)
- Gauss's averaging method for orbit-averaged torques
- Mean-motion resonance identification and libration amplitude measurement
- Kozai-Lidov cycle detection

**Planet Nine Parameters**
- Canonical parameter struct: `P9Params { mass, a, e, i, omega, Omega, M }`
- Parameter grid generation utilities
- Standard configurations from each paper era (2016 nominal, 2019 revised, 2021 MCMC)

**Initial Conditions**
- Giant planet state vectors from JPL DE440 ephemeris via `starfield`
- J2/J4 secular quadrupole field for orbit-averaged giant planet potentials
- Synthetic scattered disk generators (uniform in a, q, half-normal i)
- Inner Oort Cloud population generators

**External Perturbations**
- Galactic tidal acceleration (vertical + radial; rho_MW = 0.1 M_sun/pc^3)
- Passing star encounters (Heisler & Tremaine 1986 methodology)
- Plummer sphere cluster environment for birth cluster simulations

**I/O and Analysis**
- Orbital element time-series output (HDF5 or Parquet)
- Kernel density estimation for orbital distributions
- Phase-space portrait generation
- Snapshot sampling at configurable intervals

**Visualization**
- SVG/PNG output for orbital element distributions
- Phase-space portraits (e vs Delta-varpi)
- Mollweide sky projections for position probability
- Animated orbital evolution sequences

### Tests
- Two-body Kepler problem (verify energy conservation to machine precision)
- Circular restricted three-body problem against known Jacobi integral
- Kozai-Lidov oscillation period against analytical prediction
- Mercury6 benchmark: reproduce published test case trajectories
- Secular theory vs. direct integration comparison

---

## Crate: `p9-2016-evidence`

**Paper:** Batygin & Brown (2016) "Evidence for a Distant Giant Planet in the Solar System"
**arXiv:** 1601.05438 | **DOI:** 10.3847/0004-6256/151/2/22

### Models to Implement

**1. KBO Stability Analysis**
- 6 clones of each of 6 stable KBOs (Sedna, 2012 VP113, 2004 VN112, 2010 GB174, 2000 CR105, 2010 VZ98)
- Perturb within observational uncertainties
- Integrate 4 Gyr with four giant planets
- Output: stability classification, semi-major axis drift

**2. Analytical Secular Phase-Space Portraits**
- Secular Hamiltonian to octupole order
- Test particles at 6 semi-major axes (50–550 AU), e in (0, 0.95)
- Nominal perturber: m' = 10 M_Earth, a' = 700 AU, e' = 0.6
- Output: phase portraits of e vs. Delta-varpi showing libration islands

**3. N-body Phase-Space Portraits (Planar)**
- Parameter grid: a' in (200–2000 AU, 100 AU steps), e' in (0.1–0.9, 0.1 steps)
- Masses: 0.1, 1, 10 M_Earth
- 40 trajectories per semi-major axis bin, 4 Gyr each
- Output: anti-aligned orbit locations, mean-motion resonance identification

**4. Synthetic Scattered Disk (Planar)**
- 400 test particles, a in (50, 550) AU, q in (30, 50) AU
- Exploration suite: a' in 400–1500 AU, e' in 0.4–0.9
- 4 Gyr integration
- Output: Delta-varpi footprint, eccentricity distributions, resonance maps

**5. Synthetic Scattered Disk (3D)**
- 320 test particles, a in (150, 550) AU, half-normal i (sigma = 15 deg)
- P9: m = 10 M_Earth, a = 700 AU, e = 0.6, i = 30 deg, omega = 150 deg
- Exploratory: i' in (60–180 deg, 30 deg steps), omega' in (0–360 deg, 30 deg steps)
- 4 Gyr integration with Neptune direct + J2 for other giants
- Output: longitude of perihelion clustering, ascending node confinement, argument of perihelion vs. a

### Tests
- Reproduce Figure 1: anti-alignment of 6 observed KBOs
- Reproduce Figure 4: secular phase-space portraits at 6 semi-major axes
- Reproduce Figure 6: N-body phase portraits matching analytical theory
- Reproduce Figure 9: synthetic scattered disk with clustering
- Statistical test: P(clustering by chance) < 0.007% for 6 objects

### Visual Outputs
- Phase-space portraits (e vs. Delta-varpi) at multiple semi-major axes
- Scatter plots of final orbital elements (a vs. e, a vs. i, varpi histogram)
- Time evolution of selected test particle trajectories
- Resonance location maps

### Datasets Required
- MPC orbital elements for the 6 KBOs (publicly available from Minor Planet Center)
- JPL DE440 ephemeris for giant planet positions (via `starfield`)

---

## Crate: `p9-2016-constraints`

**Paper:** Brown & Batygin (2016) "Observational Constraints on the Orbit and Location of Planet Nine"
**arXiv:** 1603.05712 | **DOI:** 10.3847/2041-8205/824/2/L23

### Models to Implement

**1. Planar Parameter Grid (Main Survey)**
- 400 eccentric planetesimals per simulation
- a in (150, 550) AU, q in (30, 50) AU
- Parameter grid: a_9 in (200–2000 AU, 100 AU steps), e_9 in (0.1–0.9, 0.1 steps)
- Masses: 0.1, 1, 10, 20, 30 M_Earth
- Total: ~320 simulations per mass x 5 masses = 1600 simulations
- Timestep: 1/10 Neptune's orbital period
- Duration: 4 Gyr each
- Metrics: perihelion clustering probability, high-perihelion fraction

**2. 3D Inclination Exploration**
- Fixed: a_9 = 700 AU, e_9 = 0.6
- i_9 in {1, 10, 20, 30, 60, 90, 120, 150 deg}
- omega_9 in (0–360 deg, 30 deg steps) for each i_9
- 400 particles per simulation
- Metrics: confinement probability (7 random objects from a in [300,700], q < 80, i < 50, survival > 3 Gyr)
- Pole angle: average > 20 deg, RMS spread < 6.2 deg

**3. Detection Limit Assessment**
- Predicted V magnitude: 22–25 for allowed parameter space
- Map excluded sky regions from WISE, CRTS, Pan-STARRS, DES, Cassini
- Output: ~2/3 of orbit ruled out

### Tests
- Reproduce Figure 2: accepted parameter space in (a_9, e_9) plane for each mass
- Reproduce Figure 5: inclination constraints (accepted i_9, omega_9 combinations)
- Verify: mass 5–20 M_Earth, perihelion 150–350 AU, a 380–980 AU, i 22–40 deg

### Visual Outputs
- Heat maps of acceptance probability in (a_9, e_9) space per mass
- Accepted parameter volume as function of m_9
- Sky map showing excluded regions from surveys
- Corner plot of allowed (m_9, a_9, e_9, i_9) parameter ranges

### Datasets Required
- Survey footprints and depth maps for WISE, CRTS, Pan-STARRS DR1, DES Y1-Y3
- Cassini ranging residual constraints on distant massive bodies

---

## Crate: `p9-2016-obliquity`

**Paper:** Bailey, Batygin & Brown (2016) "Solar Obliquity Induced by Planet Nine"
**arXiv:** 1607.03963 | **DOI:** 10.3847/0004-6256/152/5/126

### Models to Implement

**1. Secular Hamiltonian Integration**
- Two massive wires: combined giant planet angular momentum + Planet Nine
- Solar spin vector as third body (initially aligned with giant planets)
- Poincare action-angle coordinates
- Quadrupole-order secular Hamiltonian
- Skumanich spin-down: omega_sun proportional to 1/sqrt(t)
- Solar interior: n=3 polytrope, I = 0.08, k_2 = 0.01
- Duration: 4.5 Gyr

**2. Parameter Space Exploration**
- m_9 in {10, 15, 20 M_Earth}
- a_9 in (400–900+ AU)
- e_9 in (0.3–0.9)
- Constraint: 150 < q_9 < 350 AU
- Target: produce 6 deg solar obliquity

**3. Backward Integration**
- From present obliquity, integrate backward to estimate primordial value
- Test sensitivity to initial conditions

### Tests
- Reproduce Figure 2: contours of required i_9 to achieve 6 deg obliquity in (a_9, e_9) plane
- Reproduce Figure 3: Delta-Omega contours
- Verify optimal: m_9 = 15 M_Earth, a_9 = 500 AU, e_9 = 0.5, i_9 = 20 deg
- Energy conservation check: secular Hamiltonian conserved to < 10^-10

### Visual Outputs
- Contour maps of required i_9 in (a_9, e_9) space for each mass
- Time evolution of solar obliquity over 4.5 Gyr
- Phase portraits of spin-orbit coupling
- Delta-Omega contour maps

### Datasets Required
- Solar rotation period and moment of inertia (standard solar model values)
- Giant planet orbital elements and masses (JPL DE440 via `starfield`)

---

## Crate: `p9-2016-inclined-tnos`

**Paper:** Batygin & Brown (2016) "Generation of Highly Inclined Trans-Neptunian Objects by Planet Nine"
**arXiv:** 1610.04992 | **DOI:** 10.3847/2041-8205/833/1/L3

### Models to Implement

**1. Full 4 Gyr N-body Simulation**
- 3,200 test particles
- a in (150, 550) AU, q in (30, 50) AU
- Half-normal inclination (sigma_i = 15 deg), random omega, Omega, M
- All four giant planets in direct N-body (no J2 averaging)
- P9: m = 10 M_Earth, a = 600 AU, e = 0.5, i = 30 deg, omega = 150 deg
- Timestep: 300 days
- Removal: r < 5 AU or r > 10,000 AU

**2. Dynamical Pathway Analysis**
- Track individual particle histories showing Kozai-Lidov + Neptune scattering pathway
- Identify: P9 pumps inclination via Kozai-Lidov, Neptune scattering reduces semi-major axis
- Two-step pathway produces high-i, low-a objects

### Tests
- Reproduce Figure 1: density histogram in (a, i) and (q, i) space
- Match known objects: Drac (a=41, e=0.5, i=103 deg), Niku (a=36, e=0.3, i=110 deg), 2016 NM56 (a=74, e=0.9, i=144 deg)
- Verify bimodal inclination distribution: prograde (0–110 deg) + retrograde (~150 deg)
- Reproduce Figure 2: time evolution of selected particle (a, e, i vs t)

### Visual Outputs
- 2D density histograms in (a, i) and (q, i) space
- Time-series plots of example particle trajectories showing Kozai-Lidov oscillations
- Scatter plot of surviving population colored by dynamical pathway
- Comparison overlay with observed high-inclination TNOs

### Datasets Required
- MPC orbital elements for Drac, Niku, 2016 NM56 and other high-i TNOs
- JPL DE440 ephemeris for giant planets

---

## Crate: `p9-2017-bias`

**Paper:** Brown (2017) "Observational Bias and the Clustering of Distant Eccentric Kuiper Belt Objects"
**arXiv:** 1706.04175 | **DOI:** 10.3847/1538-3881/aa79f4

### Models to Implement

**1. Observational Bias Function Computation**
- For each KBO: compute bias function B_j over (i, varpi, Omega) given (a, e, H)
- Use full MPC discovery catalog (1248 objects as of paper date)
- Smooth to 1 deg resolution
- Account for survey pointing, depth, and detection efficiency

**2. Statistical Clustering Test**
- 10,000 Monte Carlo iterations with uniform (varpi, Omega) and sin(i)*exp(-i^2/2*sigma^2) for inclinations
- Compute clustering probability for 10 KBOs with a > 230 AU
- Joint probability of longitude of perihelion AND argument of perihelion clustering

### Tests
- Reproduce clustering probability: 1.2% (varpi alone), 0.025% (combined)
- Verify bias function shapes for individual KBOs against published figures
- Sensitivity test: vary number of KBOs and a threshold

### Visual Outputs
- Bias function maps for each KBO (sky projection)
- Monte Carlo probability distribution with observed value marked
- Cumulative distribution of clustering statistic

### Datasets Required
- Full MPC discovery catalog with survey metadata (publicly available from MPC)
- Orbital elements of 10+ distant KBOs with a > 230 AU

---

## Crate: `p9-2018-kuiper-belt`

**Paper:** Khain, Batygin & Brown (2018) "The Generation of the Distant Kuiper Belt by Planet Nine from an Initially Broad Perihelion Distribution"
**arXiv:** 1804.11281 | **DOI:** 10.3847/1538-3881/aac212

### Models to Implement

**1. Narrow Initial Perihelion Distribution**
- q in (30, 36) AU, a in (150, 550) AU
- P9: m = 10 M_Earth, a = 700 AU, e = 0.6
- 4 Gyr mercury6 integration
- Output: perihelion distance distribution, alignment population ratios

**2. Broad Initial Perihelion Distribution**
- q in (30, 300) AU, a in (150, 550) AU
- Same P9 parameters
- 4 Gyr integration
- Output: bimodal perihelion structure (low-q anti-aligned + high-q aligned)

### Tests
- Verify only broad distribution produces both aligned AND anti-aligned populations
- Reproduce Figure 2: final perihelion distributions for narrow vs. broad cases
- Identify two permanently stable populations

### Visual Outputs
- Side-by-side comparison: narrow vs. broad initial q distributions
- Final (a, q) scatter plots colored by alignment/anti-alignment
- Perihelion distance histograms at 0, 1, 2, 3, 4 Gyr

### Datasets Required
- No external datasets beyond standard giant planet ephemeris

---

## Crate: `p9-2018-resonance`

**Paper:** Bailey, Brown & Batygin (2018) "Feasibility of a Resonance-based Planet Nine Search"
**arXiv:** 1809.02594 | **DOI:** 10.3847/1538-3881/aaccf4

### Models to Implement

**1. Planar Resonance Characterization**
- a_9 = 600 AU, e_9 from 0 to 0.7
- Grid of test particle semi-major axes
- Identify all mean-motion resonances (N/1, N/2, N/3, etc.)
- Measure resonance widths and libration amplitudes

**2. Resonance Probability Analysis**
- Compute probability that all 6 observed KBOs reside in N/1 or N/2 resonances
- Full resonance distribution analysis showing broad plateau vs. discrete peaks

### Tests
- Reproduce: probability < 5% that all 6 KBOs in N/1 or N/2 resonances
- Verify: a_9 ~ 660 AU peak dissolves into broad plateau with full resonance set
- Reproduce resonance width vs. eccentricity curves

### Visual Outputs
- Resonance location map in (a_particle, a_P9) space
- Libration amplitude vs. semi-major axis
- Probability density of KBO semi-major axes given resonance structure

### Datasets Required
- Orbital elements of 6 reference KBOs from MPC

---

## Crate: `p9-2019-clustering`

**Paper:** Brown & Batygin (2019) "Orbital Clustering in the Distant Solar System"
**arXiv:** 1901.07115 | **DOI:** 10.3847/1538-3881/aaf051

### Models to Implement

**1. Extended Bias Analysis**
- 14 KBOs with a > 230 AU (expanded from 10 in 2017 paper)
- Simultaneous longitude of perihelion AND orbital pole position biases
- Full MPC survey metadata

**2. OSSOS Survey Sensitivity Test**
- 4 distant KBOs from OSSOS survey
- Demonstrate that OSSOS's limited sky coverage cannot detect clustering even if real
- Survey simulator for OSSOS footprint

### Tests
- Reproduce: combined clustering probability = 0.2% (99.8% confidence)
- Verify OSSOS insensitivity claim
- Reproduce Figure 3: bias-corrected clustering statistics

### Visual Outputs
- Updated bias-corrected varpi distribution
- OSSOS footprint vs. KBO population sky map
- Probability distribution from Monte Carlo with 99.8% threshold marked

### Datasets Required
- MPC catalog (updated to 2019)
- OSSOS survey characterization files (publicly available from OSSOS team)

---

## Crate: `p9-2019-review`

**Paper:** Batygin, Adams, Brown & Becker (2019) "The Planet Nine Hypothesis"
**arXiv:** 1902.10103 | **DOI:** 10.1016/j.physrep.2019.01.009

### Models to Implement

**1. Semi-averaged Integrations (1,134 simulations)**
- Secular + octupole Hamiltonian with orbit averaging
- Large parameter grid: m_9, a_9, e_9, i_9
- Fast exploration of parameter space

**2. Fully Resolved N-body (660 simulations)**
- Direct mercury6-style integration
- Subset of parameter grid for validation
- 4 Gyr each

**3. Revised Parameter Estimation**
- From (10 M_Earth, a=700, e=0.6) to revised (5-10 M_Earth, a=400-800, e=0.2-0.5, i=15-25 deg)
- Predicted V magnitude: 19-24

### Tests
- Reproduce Table 1: parameter constraints comparison (2016 vs. 2019)
- Reproduce Figure 14: revised allowed parameter space
- Verify: lower mass + closer = brighter but cannot explain solar obliquity

### Visual Outputs
- Side-by-side parameter space comparison (2016 nominal vs. 2019 revised)
- Predicted brightness map on sky
- Comprehensive figure set from review (selected key figures)

### Datasets Required
- All datasets from previous papers (accumulated)
- Updated MPC catalog

---

## Crate: `p9-2021-oort-cloud`

**Paper:** Batygin & Brown (2021) "Injection of Inner Oort Cloud Objects Into the Distant Kuiper Belt by Planet Nine"
**arXiv:** 2104.05799 | **DOI:** 10.3847/2041-8213/abee1f

### Models to Implement

**1. Phase 1: Inner Oort Cloud Formation (Birth Cluster)**
- 10^5 test particles, 4.5–12 AU heliocentric, circular coplanar
- Plummer sphere cluster: M = 1200 M_sun, radius = 0.35 pc
- Sun at half-mass radius (~0.33 pc)
- 19 stellar species (Heisler et al. 1987 mass function)
- Velocity dispersion: 1 km/s
- Stellar encounter rate: ~50/Myr
- Jupiter (5.5 AU) + Saturn (9 AU) active
- 10 independent disks of 10^4 particles, 10 Myr each
- Timestep: 100 days, accuracy 10^-12

**2. Phase 2: P9-Driven 4.5 Gyr Evolution**
- Each of 10 IOC snapshots as initial conditions
- P9: m = 5 M_Earth, a = 500 AU, e = 0.25, i = 20 deg
- Giants via J2; Galactic tide; passing stars
- Timestep: 3,000 days
- Selection: t > 2 Gyr, 40 < q < 100 AU, i < 40 deg

**3. Scattered Disk Control**
- 1000 particles, a in (100, 800) AU, q in (30, 100) AU
- Same P9 and environment
- Baseline comparison

### Tests
- Reproduce: ~20% of IOC injected into distant Kuiper belt
- Apsidal confinement: f_varpi = 67% (IOC) vs. 88% (scattered disk)
- Reproduce Figure 2: IOC evolution snapshots at 1, 2, 5, 10 Myr
- Reproduce Figure 3: final orbital element distributions

### Visual Outputs
- IOC formation snapshots (a, e scatter plots at multiple epochs)
- Phase 2 evolution: time-series of injected fraction
- Comparison: IOC-origin vs. scattered-disk-origin orbital distributions
- Apsidal alignment histograms

### Datasets Required
- Heisler et al. (1987) stellar mass function
- Galactic tidal parameters (Nesvorny et al. 2017)
- Passing star encounter rates

---

## Crate: `p9-2021-orbit`

**Paper:** Brown & Batygin (2021) "The Orbit of Planet Nine"
**arXiv:** 2108.09868 | **DOI:** 10.3847/1538-3881/ac2056

### Models to Implement

**1. 121 N-body Simulations**
- 16,800–64,900 test particles per simulation
- a in (150–500) AU, q in (30–50) AU, i in (0–25 deg)
- All four giant planets as active perturbers
- m_9 in {3, 4, 5, 6, 7, 10, 12, 20 M_Earth}
- Manual grid of (a_9, e_9, i_9)
- Duration: 4 Gyr, analyze after 1 Gyr
- Timestep: 300 days, adaptive for encounters
- Removal: r < 4.5 AU or r > 10,000 AU

**2. Kernel Density Estimation**
- 5% smoothing for a < 230 AU, linearly rising to 30% by a = 730 AU
- Probability distributions P_jk over (i, Delta-varpi, Delta-Omega) conditioned on (a, e)

**3. Gaussian Process Emulator**
- Matern kernel (nu = 1.5) interpolating 121 simulations
- Enables continuous parameter space evaluation

**4. MCMC Posterior Inference**
- emcee (affine-invariant ensemble sampler)
- 100 chains, 20,890 samples each
- Burn-in: 260 steps; thinning: every 42 steps
- Final: 49,100 uncorrelated samples
- Two priors: uniform and cluster-scattering (Frechet)

### Tests
- Reproduce posterior: m_9 = 6.2 (+2.2/-1.3) M_Earth, a_9 = 380 (+140/-80) AU
- Reproduce: i_9 = 16 +/- 5 deg, e_9 = 0.3 (+0.1/-0.1), q_9 = 300 (+85/-60) AU
- Verify clustering significance: 99.6%
- Reproduce Figure 5: MCMC posterior corner plot
- GP emulator cross-validation: leave-one-out error < 5%

### Visual Outputs
- MCMC corner plot (m_9, a_9, e_9, i_9, omega_9, Omega_9)
- Mollweide sky projection of Planet Nine position probability
- GP emulator prediction vs. direct simulation comparison
- Brightness histogram from posterior

### Datasets Required
- Reference population of 100,000 synthetic P9 objects (archived at Caltech)
- Orbital elements of ~11 reference KBOs used in likelihood
- MPC catalog

---

## Crate: `p9-2021-ztf`

**Paper:** Brown & Batygin (2021) "A Search for Planet Nine Using the Zwicky Transient Facility Public Archive"
**arXiv:** 2110.13117 | **DOI:** 10.3847/1538-3881/ac32dd

### Models to Implement

**1. Detection Limit Calibration**
- Self-calibrate using 100,000 numbered asteroids
- HEALPix gridding (NSIDE=32, ~1.8 deg^2 cells)
- Uniform V ~ 20.5 detection limit

**2. Synthetic Population Injection**
- 100,000 synthetic P9 from BB21 MCMC posterior
- Random mean anomaly M in [0, 360 deg]
- Mass-diameter: r_9 = (m_9/3M_Earth) * R_Earth
- Albedo: 0.2–0.75

**3. Orbit Linking Algorithm**
- Holman et al. (2018) method
- Great-circle motion assumption
- Heliocentric distance sampling: uniform in 1/r (Delta = 0.0001 AU^-1)
- Angular velocity clustering: Delta-alpha_dot = 5e-5 deg/day, Delta-beta_dot = 2e-5 deg/day
- Requirement: n >= 7 detections; Bernstein RMS < 2 arcsec

### Tests
- Reproduce linking efficiency: 99.66% (56,173/56,373)
- Reproduce: 56.4% of synthetic population detected (>= 7 times)
- No false positive detections in real data

### Visual Outputs
- Detection limit map (HEALPix, Mollweide projection)
- Sky map of ruled-out vs. allowed P9 locations
- Recovery efficiency vs. apparent magnitude
- Survey coverage completeness map

### Datasets Required
- ZTF public archive (June 2018 – May 2020) via IRSA
- MPC numbered asteroid catalog for calibration
- HEALPix library

---

## Crate: `p9-2022-des`

**Paper:** Belyakov, Bernardinelli & Brown (2022) "Limits on the Detection of Planet Nine in the Dark Energy Survey"
**arXiv:** 2203.07642 | **DOI:** 10.3847/1538-3881/ac5c56

### Models to Implement

**1. DESTNOSIM Survey Simulator**
- 100,000 synthetic P9 from BB21 posterior
- DES footprint: ~5000 deg^2 southern sky
- 575 nights (2013–2019), five bands (grizY)
- Completeness: logit function p(m) = c/(1 + exp(k(m - m_50))), m_50 = 24.1

**2. Five Color/Albedo Models**
- Fiducial (BB21): g-r=0.44, g-i=0.53, g-z=0.55, albedo 0.2-0.75
- Neptune-like: albedo 0.5, g-r=-0.3, g-i=-1.35, g-z=-2.15
- 40K 0.1 CH4: albedo 0.75, g-r=0.6, g-i=-0.3, g-z=-0.2
- Super-Ganymede: albedo 0.43, g-r=0.72, g-i=0.88, g-z=0.87
- Super-KBO: albedo 0.1, g-r=1.0, g-i=1.25, g-z=1.5

**3. Recovery Analysis**
- NUNIQUE >= 7, ARCCUT >= 180 days, at least one triplet within 90 days
- H magnitude: H = m_sun - 2.5*log10(p*D^2/(9e16 km^2))

### Tests
- Reproduce: 87.0% recovery (fiducial), range 76.8–88.1% across models
- Reproduce: DES alone rules out ~10.2%; 5.0% uniquely beyond ZTF
- Combined ZTF+DES: 61.2%

### Visual Outputs
- DES footprint with P9 detection probability overlay
- Recovery fraction vs. apparent magnitude for each color model
- Venn diagram: ZTF vs. DES parameter space exclusion
- Color-dependent detection efficiency comparison

### Datasets Required
- DES survey metadata and pointing history (DES Data Management)
- DESTNOSIM survey simulator code (Bernardinelli et al. 2022)

---

## Crate: `p9-2024-panstarrs`

**Paper:** Brown, Holman & Batygin (2024) "A Pan-STARRS1 Search for Planet Nine"
**arXiv:** 2401.17977 | **DOI:** 10.3847/1538-3881/ad24e9

### Models to Implement

**1. PS1 Data Processing Pipeline**
- Query PS1 CasJobs for single-night transients
- Quality cuts: PSFQF > 0.99, PSFCHI2 in [0-8 bright, 0-1.6 faint], mag < 22.5
- Bright star masking: r_e = 50 + 0.08*(13-m_r)^4 arcsec
- Clustering removal: groups >= 5 objects within 40 arcsec
- Result: 244 million objects from 1.26 billion raw

**2. Self-Calibration**
- 501,090 numbered asteroids; 6.21M detections
- RMS offset 0.2 mag

**3. Orbit Linking**
- Holman et al. (2018) variant
- Heliocentric distance: Delta(1/r) = 10^-5 AU^-1
- Angular velocity box: 0.052 deg/day (lon), 0.026 deg/day (lat)
- Requirement: n >= 9 detections
- Processing: ~50 3x3 deg blocks in parallel

**4. Injection Test**
- 100,000 reference population
- 68,550/69,082 correctly linked (99.2%)

### Tests
- Reproduce: 69,802/100,000 detected by PS1; 17,054 unique to PS1
- Combined ZTF+DES+PS1: 78% ruled out
- Updated P9: a_9 = 500 (+170/-120) AU, V = 22.0 (+1.1/-1.4)

### Visual Outputs
- PS1 survey depth map
- Combined exclusion map (ZTF + DES + PS1)
- Updated posterior after PS1 constraints
- Remaining allowed sky region for P9

### Datasets Required
- Pan-STARRS1 DR2 via MAST/CasJobs
- MPC numbered asteroid catalog
- BB21 reference population

---

## Crate: `p9-2024-neptune-crossing`

**Paper:** Batygin, Morbidelli, Brown & Nesvorny (2024) "Generation of Low-Inclination, Neptune-Crossing TNOs by Planet Nine"
**arXiv:** 2404.11594 | **DOI:** 10.3847/2041-8213/ad3cd2

### Models to Implement

**1. P9-Inclusive 4 Gyr Integration**
- Initial conditions: Nesvorny et al. (2023) cluster_2 snapshot at t = 300 Myr
- ~2000 test particles with q > 30 AU, 100 < a < 5000 AU
- Remove all Neptune-crossing objects from ICs
- P9: m = 5 M_Earth, a = 500 AU, e = 0.25, i = 20 deg
- All four giant planets + Galactic tide + passing stars
- Timestep: 100 days, accuracy 10^-11
- Inner absorbing boundary: 1 AU; outer: 100,000 AU
- Duration: 4 Gyr
- Sample: 1 Myr intervals in final Gyr; select a > 100, q < 30, i < 40

**2. P9-Free Null Model**
- Nesvorny cluster_2 simulation extended to completion (no P9)
- Final Gyr sampled at 20 Myr intervals
- 2000-3000 orbital footprints

**3. Statistical Hypothesis Test**
- 17 well-characterized multi-opposition TNOs (a > 100, i < 40, q < 30)
- KS uniformity test on CDF values xi_j = CDF_rj(q_j)
- Logarithmic statistic zeta = product of log(xi_j)
- 10^6 Monte Carlo realizations for expected distribution

**4. OSSOS Validation**
- 100,000 synthetic TNOs: boxcar q in [15, 30], Gaussian i (sigma=15), power-law a (proportional to a^-3/2)
- OSSOS survey simulator (Lawler et al. 2018)
- 20 observations simulated

### Tests
- Reproduce: P9 model zeta = -7.9 (0.4 sigma, p = 0.41, consistent)
- Reproduce: P9-free model zeta = -16.5 (~5 sigma, p = 0.0034, rejected)
- OSSOS validation: KS p = 0.72, zeta = -9.8 (0.6 sigma)
- Reproduce Figure 1: perihelion distributions (P9 vs. P9-free)
- Reproduce Figure 2: zeta statistic distributions

### Visual Outputs
- Perihelion distribution comparison: P9 model vs. P9-free vs. observed
- CDF uniformity plots for both models
- Zeta statistic distribution with observed value marked
- (a, q, i) scatter plots of Neptune-crossing population

### Datasets Required
- Nesvorny et al. (2023) cluster_2 simulation output (request from authors or reproduce)
- 17 multi-opposition TNO orbital elements from MPC
- OSSOS survey simulator (publicly available)
- Galactic tidal parameters
