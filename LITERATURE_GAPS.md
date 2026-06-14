# Literature Gap Analysis — Papers Not Yet Ingested

Three-way sweep (June 2026): full arXiv keyword sweep 2014–2026, a targeted hunt for
critiques/alternatives/recent work, and a Semantic Scholar forward-citation crawl over the 21
implemented papers (1,016 citation edges → 454 unique citing papers, ranked by "in-degree" =
how many of our 21 they cite). The repo's 21 crates and their arXiv IDs are listed at the end.

## The headline blind spot

The codebase reproduces the **pro-P9 (Caltech) side of the debate almost exclusively** — the
clustering critiques it responds to, the alternative explanations, and the independent orbit
estimates are all absent. These are also the most reproducible papers in the field
(survey-simulator statistics and secular theory squarely within existing p9-core machinery).

## Top crate candidates (ranked by centrality × reproducibility × balance)

| Rank | arXiv | Paper | In-deg | Why |
|---|---|---|---|---|
| 1 | 2102.05601 | Napier et al. 2021 — "No Evidence for Orbital Clustering in the ETNOs" | 5 | THE published critique (joint OSSOS+DES+S&T selection functions → 17–94% consistency with uniform). Our 2019-clustering/2025-clustering crates implement the rebuttals to this paper without implementing the paper. |
| 2 | 2005.05326 + 2105.01065 | Clement & Sheppard 2020/2021 — orbital precession constraints; Neptune's distant MMRs under P9 | **10 both** | Highest citation centrality of any un-ingested work — they cite essentially our whole 2016–2019 core. Clean N-body constraint papers. |
| 3 | 1706.05348 | Shankman et al. 2017 — OSSOS VI detection biases | — | The original bias critique (predates much of our set, so low in-degree, but foundational). Survey-simulator statistics. |
| 4 | 2604.25990 | Siraj, Chyba & Tremaine 2026 — "Measuring Apsidal Clustering" | — | The live state of the art (Apr 2026): conditional-likelihood estimator; expanding the stable-ETNO sample 21→25 drops significance 2.7σ → **1.9σ**. No Caltech rebuttal in print yet. LSST-footprint-ready statistic. |
| 5 | 2410.18170 | Siraj, Chyba & Tremaine 2024 — "Orbit of a Possible Planet X" | 7 | Independent orbit inference (4.4 M⊕, 290 AU, i=6.8°) incompatible with Brown & Batygin 2021 — pairs directly against p9-2021-orbit. |
| 6 | 1804.06859 | Sefilian & Touma 2019 — self-gravitating TNO disk shepherding | 4 | The leading non-planet model; secular-analytic, very crate-friendly. Our p9-2024-oort-selfgrav implements the Caltech counter without the original. |
| 7 | 2310.20614 | Huang & Gladman 2024 — primordial sednoid alignment | 6 | Sednoid apsides converge 4.5 Gyr ago → rogue-planet-since-ejected explanation; a *present* P9 likely incompatible. Cheap secular back-integration. |
| 8 | 1602.06116 (+ INPOP19a 2020) | Fienga et al. — Cassini ranging constraints | 5 (2020) | The canonical ephemeris constraint channel; entirely absent from the repo. Tidal perturbation on Saturn–Earth range residuals. With Holman & Payne 2016 (1603.09008/1604.03180) as siblings. |
| 9 | 2003.08901 | Bernardinelli et al. 2020 — DES ETNO isotropy tests | 7 | 12 isotropy statistics on the DES sample; complements our DES exclusion crate with the statistical side. |
| 10 | 1905.09286 | Kaib et al. 2019 — OSSOS XV scattering TNOs | 8 | Observed scattering population disfavors some P9 configs; high centrality. |
| 11 | 2505.15806 + 2508.02162 | Cheng et al. 2025 (2017 OF201) + Chen et al. 2025 (Ammonite/2023 KQ14) | — | The two 2025 discoveries that stress the hypothesis: an 830 AU dwarf planet *outside* the cluster, and a fourth sednoid misaligned with the others. Both have reproducible stability/clustering analyses, and both postdate our newest data. |
| 12 | 2506.12854 | Chen (T.) et al. 2025 — AKARI all-sky far-IR search | 6 | **Directly material to p9-2025-iras-akari**: finds the Phan et al. IRAS/AKARI candidate pair inconsistent with a single moving object. Our crate reproduces the claim but not the refutation. |
| 13 | 2304.00576 (+2303.13339) | Brown & Mathur 2023 — MOND external-field effect | 5 | Galactic EFE aligns apsides toward the galactic center; cheap quadrupole secular Hamiltonian. Counter: Vokrouhlický et al. 2403.09555 (comets). |
| 14 | 1509.08920 (+ Das & Batygin 2307.00378) | Madigan & McCourt 2016 — inclination instability | 8 (Das) | The collective-gravity alternative + the Caltech suppression counter; the linear growth-rate scalings are implementable even if full self-gravitating N-body isn't. |
| 15 | 1712.04950 | Meisner et al. 2018 — WISE/NEOWISE 3π search | 6 | Same injection-recovery methodology family as our ZTF/DES/PS1 crates; would slot into the survey-exclusion table (and starfield surveylib). |
| 16 | 1612.07774 + 1603.02196 | Millholland & Laughlin 2016; Malhotra, Volk & Wang 2016 | 4 | The resonance-prediction ancestors that p9-2018-resonance rebuts — currently we implement the rebuttal without the originals. Malhotra's is nearly pure arithmetic. |
| 17 | 2508.14156 | Siraj et al. 2025 — "Planet Y" mean-plane warp | — | Distinct hypothesis (Mercury-to-Earth mass at 100–200 AU) from a 50–80 AU mean-plane warp; Laplace-plane estimation with debiasing, very crate-friendly. |
| 18 | 1608.08772 | Sheppard & Trujillo 2016 — new ETNOs corroborating clustering | — | The key independent discovery-sample paper on the pro side. |
| 19 | 2506.02144 | Holman et al. 2025 — Pan-STARRS distant-planet search, Part 1 | — | Independent of (and complementary to) our Brown et al. PS1 crate. |
| 20 | 2602.00802 | Iorio 2026 — Saturn precession vs P9/Planet X/Planet Y | — | Cheap analytic constraint crate covering all three current hypotheses at once. |

## Full category lists

### Foundational / precursor (pro-clustering side)
- 1608.08772 Sheppard & Trujillo 2016 (new ETNOs); 1606.02294 Sheppard+ 2016 (Kuiper belt edge);
  1810.00013 Sheppard+ 2018 (2015 TG387 "Goblin", in-deg 7); 1406.0715 de la Fuente Marcos 2014;
  0712.2198 Lykawka & Mukai 2008 (pre-P9 distant planet). Trujillo & Sheppard 2014 (Nature,
  2012 VP113) has **no arXiv posting** — DOI 10.1038/nature13156.

### Dynamics modeling
- 1712.06547 Hadden+ 2018 (chaotic TNO dynamics, in-deg 5); 1806.06867 Li+ 2018 (semianalytic
  secular, in-deg 6); 1706.06609 Becker+ 2017 (resonance hopping, in-deg 4); 2010.02234 Khain+
  2020 (resonance hopping II, in-deg 7); 1607.05111 Gomes+ 2016 (independent obliquity check);
  1608.01421 Lai 2016 (analytic obliquity); 2109.13307 Anderson & Kaib 2021 (detached-belt
  inclinations, in-deg 6); 2107.06296 Oldroyd & Trujillo 2021 (perihelion gap, in-deg 7);
  1808.01248 Cáceres & Gomes 2018 (low-perihelion P9, in-deg 6); 2207.11823 Lu & Laughlin 2022
  (Uranus tilt); 1605.02473 Beust 2016; 1805.05355 Becker+ 2018 (2015 BP519, in-deg 7);
  2008.11242 Köhne & Batygin 2020; 1604.06241 de la Fuente Marcos+ 2016.

### Observational searches / constraints
- Ephemeris: 1602.06116 Fienga+ 2016; Fienga+ 2020 INPOP19a (no arXiv, A&A 640 A6, in-deg 5);
  1603.09008 + 1604.03180 Holman & Payne 2016; 1512.05288 Iorio 2016; 2602.00802 Iorio 2026.
- IR/mm: 1712.04950 + 1611.00015 Meisner (WISE); 2104.10264 Naess+ 2021 (ACT); 2111.03831
  Rowan-Robinson 2022 (the IRAS claim Phan et al. follow); 2506.12854 Chen 2025 (AKARI);
  Simons Observatory forecast 2503.00636.
- Optical: 2506.02144 Holman+ 2025 (PS1 Part 1); 2504.05473 Socas-Navarro & Trujillo 2025
  (targeted parallax search); 2010.13791 Rice & Laughlin 2020 (TESS shift-stack); 1911.03676
  Payne+ 2019 (TESS yields); 2109.03758 Bernardinelli+ 2021 (DES 6-yr, in-deg 4).
- Thermal/atmosphere models feeding detectability: 1604.07424 Fortney+ 2016; 1602.07465
  Linder & Mordasini 2016; 1602.05963 Cowan+ 2016.
- Rubin/LSST status (June 2026): first images June 2025; LSST survey proper not yet started
  (Data Preview 2 due Jul–Sep 2026); **no LSST P9 papers yet** — 2604.25990 explicitly positions
  its statistic for the LSST era; 2303.02355 Schwamb+ 2023 (observing-strategy, in-deg 6).

### Statistical critiques / independent re-analyses ← biggest gap
- 2102.05601 Napier+ 2021; 1706.05348 Shankman+ 2017 (OSSOS VI); 2003.08901 Bernardinelli+ 2020;
  1610.04251 Shankman+ 2016; 1605.06575 Lawler+ 2016; 1905.09286 Kaib+ 2019 (OSSOS XV);
  1704.02444 Volk & Malhotra 2017 (warped mean plane); 2306.12847 Beaudoin+ 2023 (OSSOS XXIX);
  2402.00266 Porter+ 2024 (9:1 MMR constraint); 2202.01693 + 2106.08369 de la Fuente Marcos
  (pro-clustering, bias-independent); 2410.18170 + 2604.25990 Siraj/Chyba/Tremaine;
  1909.09673 Kavelaars+ 2019 (in-deg 4).

### Alternative explanations
- 1804.06859 Sefilian & Touma 2019 (self-gravitating disk); 1509.08920 Madigan & McCourt 2016 +
  1805.03651 + 2004.01198 + 2004.00037 + 2304.12366 (inclination-instability series) vs
  2307.00378 Das & Batygin 2023 (suppression, in-deg 8); 1909.11090 Scholtz & Unwin 2019 (PBH) +
  2005.12280 Siraj & Loeb (LSST flares) + 2004.14192 Witten 2020; 2304.00576 Brown & Mathur +
  2303.13339 Migaszewski + 1703.06682 Paučo (MOND) vs 2403.09555 Vokrouhlický+ 2024;
  2308.13765 Lykawka & Ito 2023 (Earth-mass at 250–500 AU, in-deg 4); 2508.14156 "Planet Y";
  2310.20614 + 2505.16317 Huang & Gladman (primordial alignment; flybys unlikely);
  2007.10339 Siraj & Loeb 2020 (solar binary).

### Formation / origin
- 1603.07247 Mustill+ 2016 (capture); 1602.08496 Li & Adams 2016 (cross-sections); 1603.08010
  Bromley & Kenyon 2016; 1603.08008 Kenyon & Bromley 2016 (in-situ); 1709.00418 Parker+ 2017;
  1710.08295 Eriksson+ 2018; 1609.08614 Michaely & Loeb 2016; 2505.24093 Izidoro+ 2025
  (birth-cluster trapping, in-deg 5); 2604.19413 Hochart & Portegies Zwart 2026.

### Very recent (2025–2026), not classified above
- 2501.17129 Ribeiro+ 2025 (ecliptic comets favor ~7.5 M⊕ P9); 2509.25428 Geringer-Sameth+ 2025
  (stacking methodology); 2604.12787 Harko 2026 (DM heating of P9, in-deg 8 — quantitative but
  fringe); 2606.00612 Galiazzo & Finch 2026 (listed for completeness, likely low rigor).
- IRAS/AKARI candidate status: **no refereed confirmation or refutation of Phan et al. as of
  June 2026** beyond 2506.12854's inconsistency finding; Brown has publicly noted the implied
  orbit conflicts with P9 predictions.

## Citation-crawl notes (Semantic Scholar, June 2026)

- 21/21 covered papers fetched; 1,016 citation edges; 454 unique citing papers; deduplicated
  across arXiv/published record splits by normalized title. 2404.11594's citations live on its
  published ApJL record (the arXiv record shows zero) — merged.
- In-degree leaders (cite ~half our corpus): Clement 2020 (10), Clement & Sheppard 2021 (10),
  Trujillo 2020 review chapter (10, no arXiv), Horner 2020 review (8), Kaib 2019 (8),
  Das & Batygin 2023 (8), Harko 2026 (8).
- Known noise: one S2 record ("Brito, f(R,T) cosmology", in-deg 5) is almost certainly
  reference-extraction error. Raw per-paper JSON cached in /tmp/s2cites/ (this machine).

## Covered papers (the 21 crates)

1601.05438 (evidence), 1603.05712 (constraints), 1607.03963 (obliquity), 1610.04992
(inclined-tnos), 1706.04175 (bias), 1710.01804 (dynamics), 1804.11281 (kuiper-belt),
1809.02594 (resonance — Bailey, Brown & Batygin 2018), 1901.07115 (clustering), 1902.10103
(review), 2104.05799 (oort-cloud — Batygin & Brown 2021), 2108.09868 (orbit), 2111.00305
(stability — Batygin, Mardling & Nesvorný), 2110.13117 (ztf), 2203.07642 (des), 2404.11594
(neptune-crossing), 2405.15139 (oort-selfgrav), 2401.17977 (panstarrs), 2503.07146
(2025-clustering — Pichierri & Batygin), 2504.17288 (iras-akari — Phan et al.), 2508.10119
(perturbation — Belyakov & Batygin).
