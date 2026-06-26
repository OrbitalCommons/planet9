"""Authored text for the discussion beats: per-paper hypothesis / arguments /
conclusions, and per-prereq key learnings.

Pure data (no manim import) so the assembler can read it cheaply. Keys are the
crate directory names with dashes (e.g. ``p9-2016-evidence``); prereq keys are the
scene class names (``P01Ellipse`` ...).
"""

PAPER_TEXT = {
    # ---- 2016 ----
    "p9-2016-evidence": {
        "hypothesis": "A distant, massive ninth planet gravitationally shepherds the most distant Kuiper belt objects into aligned orbits.",
        "arguments": [
            "Six stable a>250 AU KBOs share a longitude of perihelion and orbital pole far more tightly than chance.",
            "N-body integrations show an eccentric, inclined ~10 Earth-mass planet confines apsides and raises perihelia.",
            "The same planet predicts a population of high-inclination and even retrograde objects.",
        ],
        "conclusions": [
            "The joint false-alarm probability of the clustering is about 0.007% -- very unlikely to be chance.",
            "A ~10 Earth-mass planet near a ~ 700 AU, e ~ 0.6 reproduces the observed alignment.",
        ],
    },
    "p9-2016-constraints": {
        "hypothesis": "If Planet Nine is real, its sculpting leaves quantitative, testable signatures in the distant disk.",
        "arguments": [
            "Confined apsides should show a high mean resultant length in the right semimajor-axis band.",
            "P9 should over-produce detached, high-perihelion (q > 60 AU) objects relative to a planet-free disk.",
        ],
        "conclusions": [
            "Two metrics -- apsidal confinement and a detached-object surplus -- turn the hypothesis into numbers.",
            "Acceptable P9 parameters are those that reproduce both metrics simultaneously.",
        ],
    },
    "p9-2016-cassini-ranging": {
        "hypothesis": "Planet Nine perturbs Saturn enough to shift the precisely measured Earth-Saturn range.",
        "arguments": [
            "The tidal pull of P9 adds a range residual that depends on its mass, distance, and orbital phase.",
            "Cassini measured the range to ~tens of metres, leaving little room for an unmodelled perturber.",
        ],
        "conclusions": [
            "The data most favour P9 near a true anomaly of 108-129 degrees (least disturbance).",
            "Ranging both caps the mass and points toward where on its orbit P9 currently sits.",
        ],
    },
    "p9-2016-holman-payne": {
        "hypothesis": "The absence of an anomalous Earth-Saturn range signal bounds Planet Nine's mass at each distance.",
        "arguments": [
            "A too-massive or too-close P9 would imprint a residual Cassini did not detect.",
            "The tidal residual scales as mass / distance^3, so nearby solutions are tightly constrained.",
        ],
        "conclusions": [
            "Ranging sets a distance-dependent upper bound on P9's mass.",
            "Independent of Fienga et al., it agrees that closer planets must be lighter.",
        ],
    },
    "p9-2016-iorio-precession": {
        "hypothesis": "A distant planet adds anomalous perihelion precession to the known planets.",
        "arguments": [
            "P9 would torque the inner planets, especially Saturn, adding measurable precession.",
            "Planetary ephemerides show no significant excess precession.",
        ],
        "conclusions": [
            "The lack of excess precession excludes nearby and/or massive P9 solutions.",
            "Ephemeris precession is a sensitive, model-independent constraint.",
        ],
    },
    "p9-2016-obliquity": {
        "hypothesis": "An inclined Planet Nine tilts the Sun's spin axis, explaining the ~6 degree solar obliquity.",
        "arguments": [
            "P9 torques the planets' orbital plane, which over Gyr drives a precession of the solar spin axis.",
            "The induced obliquity depends on P9's mass, distance, and inclination.",
        ],
        "conclusions": [
            "A plausible P9 naturally produces ~6 degrees of solar obliquity.",
            "A long-standing puzzle becomes a prediction rather than a coincidence.",
        ],
    },
    "p9-2016-obliquity-gomes": {
        "hypothesis": "Planet Nine tilts the entire planetary plane away from the solar equator.",
        "arguments": [
            "A massive inclined perturber secularly torques the giant planets' invariable plane.",
            "The accumulated tilt reaches the observed value for reasonable P9 parameters.",
        ],
        "conclusions": [
            "Whether framed as tilting the Sun or the planets, the ~6 degree mismatch follows from P9.",
            "It is a complementary derivation of the obliquity result.",
        ],
    },
    "p9-2016-obliquity-lai": {
        "hypothesis": "The solar obliquity from Planet Nine can be derived analytically, not just numerically.",
        "arguments": [
            "A secular / resonance-overlap treatment gives obliquity as a closed function of P9's inclination.",
            "Analytic theory exposes which parameters matter and why.",
        ],
        "conclusions": [
            "The analytic theory reproduces ~6 degrees for plausible P9 inclinations.",
            "It confirms the numerical obliquity result from first principles.",
        ],
    },
    "p9-2016-secular-resonance": {
        "hypothesis": "Distant KBOs are confined by capture into a secular resonance with Planet Nine.",
        "arguments": [
            "In secular resonance the apsidal angle librates about a fixed value rather than circulating.",
            "Libration locks the alignment in place against differential precession.",
        ],
        "conclusions": [
            "A secular resonance maintains the clustering without fine-tuned initial conditions.",
            "The libration centre and amplitude are set by P9's orbit.",
        ],
    },
    "p9-2016-resonance-prediction": {
        "hypothesis": "Planet Nine traps objects in mean-motion resonances at predictable semimajor axes.",
        "arguments": [
            "Resonant capture concentrates objects near n:1 and n:2 period ratios with P9.",
            "Those locations are computable once P9's period is assumed.",
        ],
        "conclusions": [
            "TNOs should pile up at P9's resonant distances -- a falsifiable prediction.",
            "Detecting such pile-ups would strongly support the hypothesis.",
        ],
    },
    "p9-2016-sheppard-etnos": {
        "hypothesis": "Newly discovered extreme TNOs extend the clustered sample motivating Planet Nine.",
        "arguments": [
            "Dedicated deep surveys keep finding large-a objects in the same apsidal cluster.",
            "Each new object is an independent draw testing the alignment.",
        ],
        "conclusions": [
            "The growing ETNO sample has so far reinforced the clustering signal.",
            "More objects sharpen -- but also stress-test -- the hypothesis.",
        ],
    },
    "p9-2016-inclination-instability": {
        "hypothesis": "A flattened disk of eccentric orbits is self-gravitationally unstable and tilts itself, with no planet needed.",
        "arguments": [
            "Self-gravity couples the orbits' angular momenta, driving a coherent inclination growth.",
            "The instability cones the orbits to a common inclination spontaneously.",
        ],
        "conclusions": [
            "Disk self-gravity alone can raise inclinations -- a caution against over-attributing to P9.",
            "Whether it operated depends on the early disk's mass, which is uncertain.",
        ],
    },
    "p9-2016-inclined-tnos": {
        "hypothesis": "Planet Nine can pump distant objects to very high, even retrograde, inclinations.",
        "arguments": [
            "Secular and resonant interactions with an inclined P9 excite large inclinations.",
            "Such objects are hard to produce by other known mechanisms.",
        ],
        "conclusions": [
            "P9 naturally generates the puzzling high-inclination TNO population.",
            "Observed steeply-inclined objects are a predicted fingerprint.",
        ],
    },
    "p9-2016-commensurabilities": {
        "hypothesis": "Distant TNOs sit near simple period commensurabilities with Planet Nine.",
        "arguments": [
            "Resonant capture biases periods toward small-integer ratios of P9's period.",
            "An excess of near-commensurate periods would betray that capture.",
        ],
        "conclusions": [
            "A statistical excess near integer ratios is a testable resonance signature.",
            "The current sample is small, so the test is suggestive, not decisive.",
        ],
    },
    "p9-2016-wise-coadd": {
        "hypothesis": "Stacking WISE exposures (and using the W2 band) deepens the thermal search for a cold Planet Nine.",
        "arguments": [
            "Co-adding N frames improves the limiting magnitude by ~2.5*log10(sqrt(N)).",
            "W2 (4.6 um) better captures a cold body than W1.",
        ],
        "conclusions": [
            "Stacking extends WISE's reach, but only modestly for the coldest, most distant cases.",
            "All-sky thermal coverage still leaves distant solutions open.",
        ],
    },
    "p9-2016-cowan-thermal": {
        "hypothesis": "A cold Planet Nine's spectral energy distribution peaks in the far-infrared.",
        "arguments": [
            "At ~40 K the blackbody peak (Wien's law) lands near 60-90 um.",
            "Optical and near-IR bands sit far down the spectrum's blue tail.",
        ],
        "conclusions": [
            "Far-IR surveys are the natural place to look for thermal emission.",
            "The SED guides which wavebands maximise detectability.",
        ],
    },
    "p9-2016-fortney-thermal": {
        "hypothesis": "Retained formation heat makes Planet Nine warmer (and brighter in the IR) than solar equilibrium alone.",
        "arguments": [
            "A several-Earth-mass body cools slowly and keeps an internal heat reservoir.",
            "Its effective temperature exceeds the passive equilibrium value.",
        ],
        "conclusions": [
            "P9 is warmer than a passive blackbody, easing thermal detection.",
            "Brightness predictions must include internal heat, not just reflected sunlight.",
        ],
    },
    "p9-2016-linder-evolution": {
        "hypothesis": "Thermal-evolution models set Planet Nine's present-day brightness from its age and mass.",
        "arguments": [
            "Cooling tracks give luminosity and temperature as functions of mass and age.",
            "After 4.5 Gyr a few-Earth-mass body settles near ~40 K.",
        ],
        "conclusions": [
            "Today P9 is faint but warm enough to glow in the far-IR.",
            "Evolution models tie detectability to formation history.",
        ],
    },
    # ---- 2017 ----
    "p9-2017-bias": {
        "hypothesis": "Even after carefully modelling survey selection, a real clustering signal remains.",
        "arguments": [
            "Each survey's detection probability is mapped as a function of orbital angles.",
            "The debiased object excess peaks where survey sensitivity does not.",
        ],
        "conclusions": [
            "Bias alone cannot explain the observed clustering.",
            "The signal survives a careful selection-function treatment.",
        ],
    },
    "p9-2017-dynamics": {
        "hypothesis": "Planet Nine creates a stable phase-space island in which distant orbits stay anti-aligned.",
        "arguments": [
            "The secular Hamiltonian has a libration region centred on Delta-varpi = 180 degrees.",
            "Orbits inside the island circulate around the stable point rather than wandering.",
        ],
        "conclusions": [
            "Apsidal confinement is the dynamical origin of the clustering.",
            "The island's extent depends on P9's mass and eccentricity.",
        ],
    },
    "p9-2017-ossos-bias": {
        "hypothesis": "The apparent ETNO clustering may be an artefact of where and when surveys looked.",
        "arguments": [
            "Surveys detect distant objects only near perihelion and only in their pointing footprint.",
            "An intrinsically isotropic population, viewed through that window, can look clustered.",
        ],
        "conclusions": [
            "Selection effects can manufacture clustering, so the null must be modelled carefully.",
            "It is the central counter-argument the field still debates.",
        ],
    },
    "p9-2017-resonance-hopping": {
        "hypothesis": "Distant objects are not locked forever but hop between resonances with Planet Nine.",
        "arguments": [
            "An object sticks in one resonance, escapes, and is captured by a neighbour.",
            "This random walk spreads objects in semimajor axis over time.",
        ],
        "conclusions": [
            "Resonance trapping is intermittent, not permanent.",
            "Hopping shapes the long-term distribution of distant orbits.",
        ],
    },
    # ---- 2018 ----
    "p9-2018-bp519": {
        "hypothesis": "The extreme high-inclination object 2015 BP519 is the kind of body Planet Nine produces.",
        "arguments": [
            "BP519 has i ~ 54 degrees, hard to generate without a distant perturber.",
            "Its pumping history is consistent with P9-driven inclination excitation.",
        ],
        "conclusions": [
            "A real, steeply-inclined ETNO matches a P9 prediction.",
            "Such objects are difficult to explain in a planet-free solar system.",
        ],
    },
    "p9-2018-chaotic-dynamics": {
        "hypothesis": "Overlapping mean-motion resonances with Planet Nine make distant orbits chaotic.",
        "arguments": [
            "Each resonance has a finite width in semimajor axis.",
            "When neighbouring widths overlap (Chirikov's criterion), motion becomes chaotic.",
        ],
        "conclusions": [
            "Distant orbits diffuse chaotically over billions of years.",
            "Chaos sets which orbits survive and which are cleared.",
        ],
    },
    "p9-2018-kuiper-belt": {
        "hypothesis": "Planet Nine broadens the perihelion distribution of the distant scattered disk.",
        "arguments": [
            "P9 lifts a tail of objects to high perihelion, detaching them from Neptune.",
            "A planet-free disk keeps perihelia narrowly near Neptune.",
        ],
        "conclusions": [
            "P9 produces a broad, detached high-q population.",
            "The perihelion distribution is a diagnostic of P9's influence.",
        ],
    },
    "p9-2018-low-perihelion": {
        "hypothesis": "If Planet Nine itself has a low perihelion, it reshapes the distant disk differently.",
        "arguments": [
            "A more eccentric P9 brings its perihelion nearer the giant planets.",
            "That changes which orbits get scattered versus confined.",
        ],
        "conclusions": [
            "P9's own perihelion is a free parameter that alters its sculpting.",
            "Different P9 eccentricities predict different disk structures.",
        ],
    },
    "p9-2018-resonance": {
        "hypothesis": "Specific mean-motion resonances capture distant TNOs with computable probabilities.",
        "arguments": [
            "Strong low-order resonances (2:1, 3:1) capture most readily.",
            "A resonance catalogue ranks capture likelihood.",
        ],
        "conclusions": [
            "Capture favours the strong low-order resonances.",
            "A resonant sub-population would be a smoking gun for P9.",
        ],
    },
    "p9-2018-secular-dynamics": {
        "hypothesis": "Planet Nine imposes a forced eccentricity on distant orbits in the eccentricity-vector plane.",
        "arguments": [
            "In (k, h) = e(cos, sin) Delta-varpi, P9 shifts the equilibrium off-centre.",
            "Orbits circulate or librate around that forced value.",
        ],
        "conclusions": [
            "A forced eccentricity, anti-aligned with P9, organises the orbits.",
            "The secular picture explains the confinement compactly.",
        ],
    },
    "p9-2018-wise-search": {
        "hypothesis": "An all-sky WISE W1 search constrains relatively nearby or massive Planet Nine solutions.",
        "arguments": [
            "At 3.4 um a ~40 K body is deep on the Wien tail, so W1 is shallow for cold P9.",
            "WISE nonetheless covers the whole sky.",
        ],
        "conclusions": [
            "WISE rules out close/massive P9 but leaves distant cold solutions open.",
            "Whole-sky coverage trades against depth.",
        ],
    },
    # ---- 2019 ----
    "p9-2019-clustering": {
        "hypothesis": "The ETNO alignment can be quantified with circular statistics to a defensible significance.",
        "arguments": [
            "The mean resultant length R-bar measures how concentrated the angles are.",
            "Rayleigh-type tests assign a p-value under the uniform null.",
        ],
        "conclusions": [
            "Circular statistics convert 'looks aligned' into a number.",
            "The measured concentration is significant for the current sample.",
        ],
    },
    "p9-2019-ossos-scattering": {
        "hypothesis": "Planet Nine lifts scattering-disk objects above the Neptune-scattering line into detached orbits.",
        "arguments": [
            "P9 raises perihelia, decoupling objects from Neptune.",
            "This populates the stable detached region the data show.",
        ],
        "conclusions": [
            "P9 ferries objects into the detached belt.",
            "The detached population's supply constrains P9's strength.",
        ],
    },
    "p9-2019-review": {
        "hypothesis": "The collected evidence points to a smaller, closer Planet Nine than first proposed.",
        "arguments": [
            "Updated dynamical and observational constraints favour lower mass and semimajor axis.",
            "Multiple independent lines of evidence are synthesised.",
        ],
        "conclusions": [
            "The favoured planet is ~5 Earth masses near a ~ 500 AU, e ~ 0.25.",
            "By 2019 the case had broadened even as the planet shrank.",
        ],
    },
    "p9-2019-selfgrav-disk": {
        "hypothesis": "A massive self-gravitating debris disk can shepherd ETNO apsides into alignment with no planet.",
        "arguments": [
            "Collective disk gravity drives coherent apsidal precession.",
            "Enough disk mass can mimic the confinement P9 provides.",
        ],
        "conclusions": [
            "Clustering may stem from a massive disk rather than a planet.",
            "It requires substantial surviving mass in the distant disk.",
        ],
    },
    "p9-2019-tess-yield": {
        "hypothesis": "Shift-stacking TESS full-frame images can recover a fraction of Planet Nine solutions.",
        "arguments": [
            "Co-adding ~thousands of frames pushes the depth to I_C ~ 22.",
            "Detection efficiency rolls off near that limiting magnitude.",
        ],
        "conclusions": [
            "An exoplanet survey can be re-aimed at the outer solar system.",
            "TESS recovers the brighter / nearer part of P9 parameter space.",
        ],
    },
    # ---- 2020 ----
    "p9-2020-clement-precession": {
        "hypothesis": "Planet Nine sets a coupled precession pattern that keeps distant orbits aligned as they rotate.",
        "arguments": [
            "The distant objects' apsides drift together rather than independently.",
            "Coupled precession preserves the cluster over time.",
        ],
        "conclusions": [
            "Alignment is maintained dynamically, not frozen.",
            "The precession pattern is a signature of the perturber.",
        ],
    },
    "p9-2020-des-isotropy": {
        "hypothesis": "The DES distant-object sample tests whether the orbital angles are isotropic.",
        "arguments": [
            "A two-angle correlation statistic probes clustering within DES's footprint.",
            "DES's own sample appears largely isotropic.",
        ],
        "conclusions": [
            "Within its coverage, DES sees no strong clustering signal.",
            "Footprint and selection strongly shape any such test.",
        ],
    },
    "p9-2020-resonance-hopping": {
        "hypothesis": "The distant population splits into resonant, hopping, and scattering classes.",
        "arguments": [
            "Classification by dynamical behaviour quantifies each class's fraction.",
            "A resonant sub-population would point to P9.",
        ],
        "conclusions": [
            "The class fractions test the resonant fingerprint of P9.",
            "Distinguishing the classes requires good orbits and long integrations.",
        ],
    },
    "p9-2020-secular-octupole": {
        "hypothesis": "Octupole-order terms in Planet Nine's secular potential reshape the apsidal libration.",
        "arguments": [
            "Going beyond the quadrupole modulates the libration amplitude.",
            "The octupole widens the range of confined orbits.",
        ],
        "conclusions": [
            "Higher-order terms enlarge the confined zone.",
            "Accurate secular models need the octupole.",
        ],
    },
    "p9-2020-tess-shiftstack": {
        "hypothesis": "Targeted shift-stacking of TESS recovers known distant TNOs and sets P9's reach.",
        "arguments": [
            "Co-adding along the correct trial track stacks the moving signal.",
            "The wrong track smears the flux, suppressing false positives.",
        ],
        "conclusions": [
            "The method recovers real movers and bounds P9 detectability.",
            "Reach depends on the depth gain from frame stacking.",
        ],
    },
    # ---- 2021 ----
    "p9-2021-act-mm": {
        "hypothesis": "The ACT millimetre maps can be searched for a moving Planet Nine point source.",
        "arguments": [
            "A cold body emits near its blackbody peak in the mm.",
            "Non-detection in ACT constrains warm/near solutions.",
        ],
        "conclusions": [
            "CMB telescopes double as outer-solar-system surveys.",
            "ACT's sensitivity reaches several hundred AU for plausible masses.",
        ],
    },
    "p9-2021-des-catalog": {
        "hypothesis": "A well-characterised DES catalogue of distant objects underpins the clustering and exclusion tests.",
        "arguments": [
            "DES yields deep multi-band detections with known selection.",
            "The catalogue supplies orbits for statistical tests.",
        ],
        "conclusions": [
            "A characterised catalogue is the raw material for rigorous tests.",
            "Its depth and footprint set what can be concluded.",
        ],
    },
    "p9-2021-detached-inclinations": {
        "hypothesis": "Planet Nine forces broader, higher inclinations in the detached Kuiper belt.",
        "arguments": [
            "P9 pumps the orbital planes of detached objects.",
            "A planet-free model predicts a narrower inclination spread.",
        ],
        "conclusions": [
            "The detached belt's inclinations favour a P9-forced model.",
            "Inclination structure is an additional P9 diagnostic.",
        ],
    },
    "p9-2021-napier-critique": {
        "hypothesis": "Accounting for survey selection, the ETNO data do not require Planet Nine.",
        "arguments": [
            "Modelling the pointed surveys' biases reproduces the apparent clustering.",
            "The significance drops once selection is included.",
        ],
        "conclusions": [
            "The clustering may be consistent with a uniform parent population.",
            "Careful, survey-by-survey debiasing is essential before claiming a signal.",
        ],
    },
    "p9-2021-oort-cloud": {
        "hypothesis": "Planet Nine may be an ice giant scattered into the inner Oort cloud during the Sun's birth cluster.",
        "arguments": [
            "Stellar encounters in a dense cluster can park bodies at a ~ 2000-20000 AU.",
            "The same process supplies detached objects.",
        ],
        "conclusions": [
            "P9's origin may lie in the Sun's birth environment.",
            "The inner Oort cloud links P9 to the detached population.",
        ],
    },
    "p9-2021-orbit": {
        "hypothesis": "A Markov-chain fit to all constraints yields a posterior on Planet Nine's orbit.",
        "arguments": [
            "Dynamical and survey constraints are combined statistically.",
            "The posterior median lands near 6 Earth masses, a ~ 380 AU, e ~ 0.3.",
        ],
        "conclusions": [
            "The posterior narrows where to point next.",
            "It seeds this project's reference population for every survey-reach figure.",
        ],
    },
    "p9-2021-perihelion-gap": {
        "hypothesis": "A deficit of objects at intermediate perihelion separates the scattering and detached populations.",
        "arguments": [
            "Objects are scarce at q ~ 42-55 AU between the two regimes.",
            "Planet Nine can carve such a gap.",
        ],
        "conclusions": [
            "The perihelion gap is a structural fingerprint of the distant disk.",
            "Its presence is consistent with P9-driven sculpting.",
        ],
    },
    "p9-2021-stability": {
        "hypothesis": "Distant scattered-disk orbits survive 4.5 Gyr only outside a perihelion/semimajor-axis boundary.",
        "arguments": [
            "Inside the boundary, objects are cleared by repeated encounters.",
            "Only sufficiently detached orbits persist.",
        ],
        "conclusions": [
            "The stability boundary bounds the dynamically allowed P9 family.",
            "Survival over the age of the system is a hard constraint.",
        ],
    },
    "p9-2021-ztf": {
        "hypothesis": "A direct ZTF shift-and-link search would have detected the brightest Planet Nine solutions.",
        "arguments": [
            "ZTF reaches r ~ 20.5 over most of the northern sky.",
            "No moving source consistent with P9 was found.",
        ],
        "conclusions": [
            "ZTF rules out about 56.4% of the predicted parameter space.",
            "More than half the hiding places are now eliminated.",
        ],
    },
    # ---- 2022 ----
    "p9-2022-des": {
        "hypothesis": "The deep but narrow DES survey rules out faint, distant solutions in its footprint.",
        "arguments": [
            "DES reaches r ~ 23.8 but covers only ~5000 deg^2 (~12% of sky).",
            "Exclusion only applies to objects that fall in that footprint.",
        ],
        "conclusions": [
            "DES adds ~5% unique exclusion after ZTF.",
            "Depth buys little if the footprint is small -- coverage is king.",
        ],
    },
    "p9-2022-iras-candidate": {
        "hypothesis": "A moving far-IR source in the IRAS archive could be a cold, distant Planet Nine.",
        "arguments": [
            "IRAS 60 um flux consistent with a cold body at a few hundred AU is the selection.",
            "A single epoch flags candidates but cannot confirm motion.",
        ],
        "conclusions": [
            "The archive yields tantalising single-epoch flags.",
            "A second epoch is required to confirm proper motion.",
        ],
    },
    "p9-2022-uranus-tilt": {
        "hypothesis": "A primordial spin-orbit resonance sweep can tip a planet's spin axis far over.",
        "arguments": [
            "A slow resonance crossing in the early solar system drives large obliquity change.",
            "It is the same obliquity machinery invoked for the Sun under P9.",
        ],
        "conclusions": [
            "Resonant obliquity sweeps are a real, demonstrated mechanism.",
            "They lend plausibility to the P9 solar-tilt scenario.",
        ],
    },
    # ---- 2023 ----
    "p9-2023-lsst-strategy": {
        "hypothesis": "Rubin/LSST's discovery fraction for Planet Nine is set by its observing-strategy choices.",
        "arguments": [
            "Per-visit depth, nights-to-link, and footprint reach each drive the yield.",
            "Tuning these levers changes how much of the viable region is searched.",
        ],
        "conclusions": [
            "A well-designed LSST strategy discovers most of the viable region.",
            "Survey design, not just telescope size, decides the outcome.",
        ],
    },
    "p9-2023-mond-efe": {
        "hypothesis": "Modified Newtonian Dynamics' external-field effect can align distant apsides without a planet.",
        "arguments": [
            "The galaxy's external field imposes a weak, fixed-direction tide on very distant orbits.",
            "Over Gyr this aligns apsides, producing clustering.",
        ],
        "conclusions": [
            "Clustering might require modified gravity rather than a ninth planet.",
            "It is a serious alternative to the P9 hypothesis.",
        ],
    },
    # ---- 2024 ----
    "p9-2024-neptune-crossing": {
        "hypothesis": "Some distant objects dip to Neptune-crossing perihelia, testing the P9 and scattering picture.",
        "arguments": [
            "Neptune-crossers are short-lived and must be resupplied.",
            "Their numbers constrain the distant disk's structure.",
        ],
        "conclusions": [
            "The supply rate of crossers bounds the distant population.",
            "It connects the observed TNOs to their dynamical source.",
        ],
    },
    "p9-2024-oort-selfgrav": {
        "hypothesis": "The inner Oort cloud's own self-gravity drives coherent apsidal precession.",
        "arguments": [
            "Collective self-gravity couples the orbits' precession.",
            "It shapes the detached population independent of, or with, P9.",
        ],
        "conclusions": [
            "A massive cloud can stir itself into coherent motion.",
            "Self-gravity is another route to apsidal organisation.",
        ],
    },
    "p9-2024-panstarrs": {
        "hypothesis": "Adding a Pan-STARRS1 search to ZTF and DES sharply shrinks the viable Planet Nine space.",
        "arguments": [
            "PS1 3-pi covers the sky north of dec -30 to r ~ 21.5 over many epochs.",
            "The union of independent searches is taken per orbit.",
        ],
        "conclusions": [
            "The combined searches rule out ~78.5% of the prior.",
            "About one-fifth of the parameter space survives -- where Rubin will look.",
        ],
    },
    "p9-2024-primordial-alignment": {
        "hypothesis": "The sednoids' apsides may have been aligned at birth and simply preserved.",
        "arguments": [
            "A stellar flyby or the birth cluster could imprint an initial alignment.",
            "Rewinding the orbits converges to a tighter primordial cluster.",
        ],
        "conclusions": [
            "Alignment might be a frozen relic, needing no present-day planet.",
            "It is an alternative origin for the clustering.",
        ],
    },
    "p9-2024-siraj-orbit": {
        "hypothesis": "Combining all constraints yields a refined posterior on Planet Nine's orbit.",
        "arguments": [
            "Dynamical and survey limits are folded into a single estimate.",
            "The result narrows mass, semimajor axis, eccentricity, and inclination.",
        ],
        "conclusions": [
            "The best estimate is a few-Earth-mass body near ~450 AU.",
            "A tighter posterior focuses the observational search.",
        ],
    },
    # ---- 2025 ----
    "p9-2025-akari-refutation": {
        "hypothesis": "The proposed IRAS/AKARI far-IR candidate is not actually Planet Nine.",
        "arguments": [
            "The candidate's implied orbit conflicts with dynamical constraints.",
            "Its motion/flux does not hold up under scrutiny.",
        ],
        "conclusions": [
            "The candidate does not survive follow-up.",
            "Extraordinary claims require corroborating evidence.",
        ],
    },
    "p9-2025-clustering": {
        "hypothesis": "Both the degree of clustering and its rate of diffusion can be measured.",
        "arguments": [
            "A statistic tracks how concentrated the angles are.",
            "Clustering diffuses with time under perturbations.",
        ],
        "conclusions": [
            "Seeing clustering today implies a force maintaining it.",
            "Diffusion sets how long alignment can persist unaided.",
        ],
    },
    "p9-2025-iras-akari": {
        "hypothesis": "A cold Planet Nine would have moved measurably between the IRAS (1983) and AKARI (2006) surveys.",
        "arguments": [
            "Cross-matching far-IR sources with consistent fluxes and the right ~23-year motion selects candidates.",
            "The implied distance is ~500-700 AU.",
        ],
        "conclusions": [
            "The pair search yields a short candidate list (one good candidate).",
            "Thermal surveys give an independent detection channel needing follow-up.",
        ],
    },
    "p9-2025-new-discoveries": {
        "hypothesis": "New distant-TNO discoveries sharpen the test of the clustering hypothesis.",
        "arguments": [
            "Objects like 2017 OF201 and Ammonite extend the sample.",
            "Some fit the P9 picture; some sit awkwardly.",
        ],
        "conclusions": [
            "Each new world is a fresh, unforgiving test.",
            "The hypothesis must accommodate the outliers as well as the fits.",
        ],
    },
    "p9-2025-parallax-search": {
        "hypothesis": "Earth's orbital baseline lets a targeted parallax search isolate a nearby body.",
        "arguments": [
            "A nearby object shows a large reflex parallax over six months.",
            "Image pairs months apart isolate that motion without long arcs.",
        ],
        "conclusions": [
            "Two epochs plus Earth's baseline can pin a nearby planet's distance.",
            "Parallax is a direct, geometric distance handle.",
        ],
    },
    "p9-2025-perturbation": {
        "hypothesis": "Perturbation theory can describe scattered-disk stability under Planet Nine analytically.",
        "arguments": [
            "An analytic series captures how P9 stabilises or destabilises orbits.",
            "It complements expensive numerical stability maps.",
        ],
        "conclusions": [
            "The analytic boundary reproduces the N-body results.",
            "Theory predicts the stability edge without giant simulations.",
        ],
    },
    "p9-2025-planet-y": {
        "hypothesis": "A separate, smaller, closer perturber ('Planet Y') would warp the Kuiper belt's mean plane.",
        "arguments": [
            "A Mercury-to-Earth-mass inclined body tilts the local Laplace plane.",
            "This signal is distinct from the distant Planet Nine.",
        ],
        "conclusions": [
            "Not every anomaly needs Planet Nine.",
            "A nearer small world can warp the belt's plane instead.",
        ],
    },
    "p9-2025-ps1-holman": {
        "hypothesis": "A pixel-level shift-and-stack of Pan-STARRS reaches deeper than catalogue linking.",
        "arguments": [
            "Co-adding pixels along millions of trial orbits gains ~1 magnitude of depth.",
            "It probes fainter movers the catalogue search misses.",
        ],
        "conclusions": [
            "Digital tracking digs below the single-image limit.",
            "The gain comes at the cost of an enormous trial space.",
        ],
    },
    "p9-2025-simons-forecast": {
        "hypothesis": "Wide millimetre surveys (ACT, Simons Observatory) reach far beyond optical for a cold Planet Nine.",
        "arguments": [
            "At the body's mm blackbody peak, sensitivity is high.",
            "SO's deep maps push the reach to large distances.",
        ],
        "conclusions": [
            "Simons Observatory could reach ~900 AU for a 5 Earth-mass world.",
            "Millimetre surveys complement optical depth.",
        ],
    },
    "p9-2025-stacking": {
        "hypothesis": "Matched-filter digital tracking over many years of images can find a faint mover.",
        "arguments": [
            "A metric over orbit-parameter space localises a candidate peak.",
            "Vast trial counts demand an honest look-elsewhere penalty.",
        ],
        "conclusions": [
            "Stacking years of data is powerful if trials are counted honestly.",
            "Significance must be corrected for the huge search space.",
        ],
    },
    # ---- 2026 ----
    "p9-2026-apsidal-clustering": {
        "hypothesis": "A calibrated statistic can rigorously measure the ETNO apsidal clustering.",
        "arguments": [
            "An estimator quantifies the concentration of perihelion directions.",
            "Its null distribution gives a defensible significance.",
        ],
        "conclusions": [
            "The estimator turns 'they look aligned' into a number with error bars.",
            "It puts the clustering claim on firmer statistical ground.",
        ],
    },
    "p9-2026-iorio-precession": {
        "hypothesis": "An extended (distended) Planet Nine perturbs Saturn's perihelion precession differently than a point mass.",
        "arguments": [
            "A non-point-mass P9 changes the precession it induces.",
            "Ephemeris precision tightens the bound either way.",
        ],
        "conclusions": [
            "Even P9's mass distribution leaves a tiny fingerprint on Saturn.",
            "Precession remains a sensitive, model-dependent probe.",
        ],
    },
    "p9-2025-russell-albedo": {
        "hypothesis": "If Planet Nine is a ~6.6 Earth-mass mini-Neptune, its size, colour, and brightness follow from exoplanet mass-radius models.",
        "arguments": [
            "Confirmed cold exoplanets (Teq < 600 K) fix the most likely radius at 2.0-2.6 Earth radii with a thin H-He envelope.",
            "Giant-planet albedos extrapolate to a V-band geometric albedo of 0.33-0.47 for such a world.",
        ],
        "conclusions": [
            "Planet Nine should sit near absolute magnitude -6.1 to -5.2 and apparent +21.9 to +22.7.",
            "That is within reach of upcoming surveys -- and would even show a resolvable disk to the largest telescopes.",
        ],
    },
    "p9-2025-stellar-flybys": {
        "hypothesis": "A close stellar flyby early in the Solar System's life could have detached the sednoids -- no present-day planet needed.",
        "arguments": [
            "The detached extreme TNOs share a low-inclination (i < 30 degrees) profile and a primordial apsidal alignment.",
            "Only special flyby geometries (coplanar or ecliptic-symmetric) reproduce that, and close encounters (q* < 1000 AU) are rare.",
        ],
        "conclusions": [
            "Folding in how often such an encounter actually happens leaves a probability of only about 5%.",
            "An early flyby is an unlikely sole sculptor of the detached belt.",
        ],
    },
    "p9-2026-stellar-companions": {
        "hypothesis": "The Sun's nearest companion might be a hidden, gravity-only object at hundreds-to-thousands of AU rather than a star at a parsec.",
        "arguments": [
            "A tidal mass-distance envelope calibrated to Planet-Nine-like constraints grows as distance cubed.",
            "It allows Earth-to-sub-Saturn masses at 300-1000 AU, rising to a Jupiter mass near 2000 AU.",
        ],
        "conclusions": [
            "A smooth dark-matter halo cannot supply it -- only a sub-Pluto mass lies within 1000 AU.",
            "Any such companion must be a structured, bound object, visible only through its gravity.",
        ],
    },
    "p9-2026-alpha-slope": {
        "hypothesis": "How a body's brightness fades with distance (flux ~ d^alpha) reveals whether it shines by reflection or its own light.",
        "arguments": [
            "Reflected sunlight falls as d^-4; self-luminous emission falls as d^-2.",
            "Of 8,557 archival TNO photometry bins, only 186 survive the eligibility cuts -- 53 reflected, 24 self-luminous, 109 anomalous.",
        ],
        "conclusions": [
            "All 24 'self-luminous' detections trace to a single instrument -- calibration systematics, not new physics.",
            "Archival photometry is too messy to run the test cleanly even on Pluto.",
        ],
    },
}


# Prereq key learnings, keyed by scene class name.
PREFACE_LEARNING = {
    "P01Ellipse": "A bound orbit is an ellipse with the Sun at one focus; its size and shape are fixed by the semi-major axis a and eccentricity e.",
    "P02KeplerLaws": "Bodies move fastest at perihelion and slowest at aphelion, and farther orbits take longer (P^2 = a^3). Planet Nine's orbit takes ~10,000 years -- we cannot wait one.",
    "P03Elements": "Six numbers orient an orbit. The longitude of perihelion (varpi = node + argument) is the direction the orbit points -- and that is what 'clustering' is about.",
    "P04Geography": "Beyond Neptune (30 AU) lie the Kuiper belt, the scattered disk, and the extreme TNOs at hundreds of AU -- a vast, barely-mapped realm where Planet Nine would hide.",
    "P05Clustering": "The extreme TNOs point the same way far more than chance allows (a high mean resultant length). Something appears to be herding their orbits.",
    "P06Precession": "Real orbits precess: their orientations slowly drift. Because that drift differs from orbit to orbit, maintaining alignment requires an active, organising force.",
    "P07Resonance": "When orbital periods form a simple integer ratio, gravitational kicks repeat and accumulate -- letting a planet trap or protect distant orbits.",
    "P08Sculpting": "A distant, massive planet can confine apsides, lift perihelia, and tilt inclinations all at once -- explaining several odd features of the outer solar system with one cause.",
    "P09ReflectedLight": "Reflected brightness falls roughly as one over distance to the fourth power, so at ~600 AU Planet Nine is only about magnitude 22 -- faint, but reachable by deep surveys.",
    "P10Thermal": "Even at ~40 K, Planet Nine glows by its own heat, peaking in the far-infrared and millimetre. That is a second, independent way to find it.",
    "P11Surveys": "A survey only finds what its depth and footprint allow, and movers are found across epochs. Where you can look shapes what you find -- coverage can fake or hide clustering.",
    "P12Indirect": "Even unseen, a planet tugs on spacecraft, precesses the planets, and can tilt the Sun -- so dynamics constrain Planet Nine without a direct image.",
    "P13Map": "Near-all-sky surveys rule out Planet Nine below ~600 AU; deeper surveys cover little sky. A viable region survives -- and Rubin/LSST will sweep most of it soon.",
}
