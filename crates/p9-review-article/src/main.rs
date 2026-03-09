use std::fmt::Write;
use std::fs;

mod diagrams;

fn main() {
    let html = generate_review_article();
    fs::write("review_article.html", &html).expect("Failed to write review article");
    println!("Generated review_article.html ({} bytes)", html.len());
}

fn generate_review_article() -> String {
    let mut html = String::with_capacity(128_000);

    write_header(&mut html);
    write_abstract(&mut html);
    write_section_1_introduction(&mut html);
    write_section_2_initial_evidence(&mut html);
    write_section_3_constraints(&mut html);
    write_section_4_dynamical_consequences(&mut html);
    write_section_5_statistical_evidence(&mut html);
    write_section_6_parameter_evolution(&mut html);
    write_section_7_observational_campaigns(&mut html);
    write_section_8_latest_evidence(&mut html);
    write_section_9_current_state(&mut html);
    write_section_10_future(&mut html);
    write_references(&mut html);
    write_footer(&mut html);

    html
}

fn write_header(html: &mut String) {
    html.push_str(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Planet Nine: A Computational Review of the Hypothesis (2016–2024)</title>
<style>
  :root {
    --bg: #fafafa;
    --fg: #1a1a1a;
    --accent: #2c5f8a;
    --muted: #666;
    --border: #ddd;
    --code-bg: #f0f0f0;
  }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: 'Georgia', 'Times New Roman', serif;
    line-height: 1.7;
    color: var(--fg);
    background: var(--bg);
    max-width: 900px;
    margin: 0 auto;
    padding: 2rem 1.5rem;
  }
  h1 {
    font-size: 1.8rem;
    text-align: center;
    margin-bottom: 0.3rem;
    color: var(--accent);
  }
  .subtitle {
    text-align: center;
    color: var(--muted);
    font-style: italic;
    margin-bottom: 2rem;
    font-size: 1.05rem;
  }
  h2 {
    font-size: 1.35rem;
    color: var(--accent);
    margin-top: 2.5rem;
    margin-bottom: 0.8rem;
    border-bottom: 2px solid var(--border);
    padding-bottom: 0.3rem;
  }
  h3 {
    font-size: 1.1rem;
    color: var(--fg);
    margin-top: 1.5rem;
    margin-bottom: 0.5rem;
  }
  p { margin-bottom: 1rem; text-align: justify; }
  .abstract {
    background: #f5f5f5;
    border-left: 4px solid var(--accent);
    padding: 1rem 1.5rem;
    margin: 1.5rem 0;
    font-style: italic;
  }
  .figure {
    margin: 2rem 0;
    text-align: center;
  }
  .figure svg {
    max-width: 100%;
    height: auto;
  }
  .figure-caption {
    font-size: 0.9rem;
    color: var(--muted);
    margin-top: 0.5rem;
    font-style: italic;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0;
    font-size: 0.95rem;
  }
  th, td {
    border: 1px solid var(--border);
    padding: 0.5rem 0.8rem;
    text-align: center;
  }
  th {
    background: var(--accent);
    color: white;
    font-weight: 600;
  }
  tr:nth-child(even) { background: #f5f5f5; }
  .highlight {
    background: #e8f0fe;
    border: 1px solid #b8d4f0;
    padding: 1rem;
    border-radius: 4px;
    margin: 1rem 0;
  }
  .ref-list {
    font-size: 0.9rem;
    line-height: 1.5;
  }
  .ref-list p {
    text-indent: -1.5rem;
    padding-left: 1.5rem;
    margin-bottom: 0.4rem;
  }
  sup { font-size: 0.75rem; }
</style>
</head>
<body>
"#,
    );
}

fn write_abstract(html: &mut String) {
    html.push_str(
        r#"
<h1>Planet Nine: A Computational Review of the Hypothesis</h1>
<p class="subtitle">Synthesis of Numerical Models from Batygin, Brown et al. (2016–2024)</p>

<div class="abstract">
<strong>Abstract.</strong> We present a computational re-implementation and synthesis of all major
numerical models from the Planet Nine literature spanning 2016–2024. Through faithful
reproduction of 15 peer-reviewed studies in a unified Rust codebase, we trace the evolution
of the Planet Nine hypothesis from its initial proposal (m&#8776;10 M<sub>&#8853;</sub>,
a&#8776;700 AU) through successive refinements to current best estimates
(m&#8776;6.2<sup>+2.2</sup><sub>&#8722;1.3</sub> M<sub>&#8853;</sub>,
a&#8776;380<sup>+140</sup><sub>&#8722;80</sub> AU). We document how observational surveys
(ZTF, DES, Pan-STARRS1) have progressively excluded 78% of the viable parameter space,
while dynamical evidence has strengthened — most recently with the ~5&#963; rejection of
a Planet Nine-free model from Neptune-crossing TNO analysis. The remaining 22% of parameter
space, concentrated at V&#8819;22 mag, awaits exploration by the Vera C. Rubin Observatory.
</div>
"#,
    );
}

fn write_section_1_introduction(html: &mut String) {
    html.push_str(
        r#"
<h2>1. Introduction</h2>

<p>The outer solar system beyond Neptune harbors a population of small bodies whose
orbital architecture has proven difficult to explain through known gravitational
influences alone. Beginning with the discovery of Sedna in 2003 — an object with a
perihelion distance of 76 AU, well beyond Neptune's gravitational reach — evidence
has accumulated that an unseen massive body may sculpt the distant Kuiper Belt.</p>

<p>In 2016, Batygin &amp; Brown presented the first comprehensive dynamical argument
for a distant giant planet, dubbed "Planet Nine" (hereafter P9), based on the
orbital clustering of six trans-Neptunian objects (TNOs) with semi-major axes
exceeding 250 AU. Their proposal triggered an intensive multi-year campaign of
theoretical modeling, statistical analysis, and observational searches that
continues to the present day.</p>

<p>This review synthesizes findings from 15 papers spanning 2016–2024, each of
which has been computationally re-implemented as part of a unified numerical
framework. We trace the evolution of the hypothesis through four phases:
initial evidence and parameter estimation (§2–3), exploration of dynamical
consequences (§4), statistical refinement (§5–6), and observational search
campaigns (§7–8). We conclude with an assessment of the current state of the
hypothesis and prospects for definitive detection (§9–10).</p>
"#,
    );
}

fn write_section_2_initial_evidence(html: &mut String) {
    html.push_str(
        r#"
<h2>2. Initial Evidence: Orbital Clustering (2016)</h2>

<h3>2.1 The Six-Object Anomaly</h3>

<p>Batygin &amp; Brown (2016a) identified six dynamically stable Kuiper Belt objects
— Sedna, 2012 VP<sub>113</sub>, 2004 VN<sub>112</sub>, 2010 GB<sub>174</sub>,
2000 CR<sub>105</sub>, and 2010 VZ<sub>98</sub> — whose orbits share a remarkable
property: their longitudes of perihelion (&#982;) cluster around
71° &pm; 16°, their longitudes of ascending node (&#937;) cluster around
113° &pm; 13°, and their arguments of perihelion (&#969;) cluster around
318° &pm; 8°. The probability of this alignment arising by chance from a uniform
distribution is approximately 0.007% — a compelling anomaly demanding explanation.</p>

<h3>2.2 The Secular Mechanism</h3>

<p>The proposed mechanism operates through secular gravitational coupling. A distant
eccentric planet, anti-aligned with the clustered KBOs (&#916;&#982; &#8776; 180°),
creates a quadrupole potential that confines test particle perihelia. The secular
Hamiltonian, expanded to octupole order, produces phase-space portraits with
stable islands at &#916;&#982; &#8776; &#960; — particles librate around the
anti-aligned configuration rather than precessing freely.</p>

<p>N-body simulations with the four giant planets plus P9 confirmed this picture:
starting from a uniform scattered disk, the secular interaction with P9 naturally
produces the observed clustering on gigayear timescales. The nominal parameters
proposed were:</p>

<table>
<tr><th>Parameter</th><th>2016 Nominal Value</th></tr>
<tr><td>Mass (m<sub>9</sub>)</td><td>10 M<sub>&#8853;</sub></td></tr>
<tr><td>Semi-major axis (a<sub>9</sub>)</td><td>700 AU</td></tr>
<tr><td>Eccentricity (e<sub>9</sub>)</td><td>0.6</td></tr>
<tr><td>Inclination (i<sub>9</sub>)</td><td>30°</td></tr>
<tr><td>Perihelion (q<sub>9</sub>)</td><td>280 AU</td></tr>
</table>
"#,
    );
}

fn write_section_3_constraints(html: &mut String) {
    html.push_str(r#"
<h2>3. Observational Constraints (2016)</h2>

<p>Brown &amp; Batygin (2016) systematically explored the parameter space through
a 19 &times; 9 &times; 5 grid spanning semi-major axis (200–2000 AU), eccentricity
(0.1–0.9), and mass (0.1–30 M<sub>&#8853;</sub>). Each point was evaluated against
three acceptance criteria: (1) at least 7 surviving test particles in the
a &#8712; [300, 700] AU range, (2) clustering fraction exceeding 50%, and
(3) generation of high-perihelion objects (q &gt; 60 AU).</p>

<p>The viable region forms a diagonal band in (a<sub>9</sub>, e<sub>9</sub>) space,
reflecting the constraint that the perihelion distance q<sub>9</sub> = a<sub>9</sub>(1 &#8722; e<sub>9</sub>)
must fall within ~200–400 AU for effective secular coupling. Too close and P9
ejects the KBOs; too far and the coupling timescale exceeds the age of the
solar system.</p>

<p>Brightness estimates using the mass-radius relation for sub-Neptune planets
(R &#8733; M<sup>0.55</sup>) and geometric albedo p &#8776; 0.3–0.75 predicted
apparent V-band magnitudes of 22–25, explaining why P9 has eluded detection by
existing all-sky surveys with limiting magnitudes of ~20–21.</p>
"#);
}

fn write_section_4_dynamical_consequences(html: &mut String) {
    html.push_str(
        r#"
<h2>4. Dynamical Consequences</h2>

<h3>4.1 Solar Obliquity</h3>

<p>Bailey, Batygin &amp; Brown (2016) demonstrated that Planet Nine provides a
natural explanation for the 6° misalignment between the Sun's spin axis and the
invariable plane of the solar system. Through secular spin-orbit coupling between
P9's inclined orbit, the giant planet orbital plane, and the solar spin axis
(modeled as an n = 3 polytrope with dimensionless moment of inertia &#206; = 0.08
and Love number k<sub>2</sub> = 0.01), a P9 inclination of 15–25° produces the
observed obliquity over the age of the solar system.</p>

<p>The mechanism requires accounting for the Sun's spin-down via Skumanich's law
(&#969;<sub>&#8857;</sub> &#8733; t<sup>&#8722;1/2</sup>), which reduces the solar
quadrupole moment and spin angular momentum over time, allowing the giant planet
plane to precess away from the solar equator. This represents a prediction unique
to the P9 hypothesis — no other proposed explanation simultaneously accounts for
both the KBO clustering and the solar obliquity.</p>

<h3>4.2 Highly Inclined and Retrograde TNOs</h3>

<p>Batygin &amp; Brown (2016b) showed that P9 also explains the existence of
highly inclined (i &gt; 60°) and retrograde (i &gt; 90°) trans-Neptunian objects
such as Drac (2008 KV<sub>42</sub>, i &#8776; 103°), Niku (2011 KT<sub>19</sub>,
i &#8776; 110°), and 2016 NM<sub>56</sub> (i &#8776; 144°). These objects are
produced through a two-step process:</p>

<p>(1) Kozai-Lidov oscillations driven by P9 pump the inclination of distant
scattered disk objects to extreme values, exchanging eccentricity for inclination
while conserving the Kozai integral &#8730;(1 &#8722; e<sup>2</sup>) cos i.</p>

<p>(2) Neptune scattering subsequently reduces the semi-major axis, decoupling the
objects from P9's influence and stranding them on stable, high-inclination orbits
in the a &#8776; 30–80 AU range.</p>

<h3>4.3 Kuiper Belt Bimodality</h3>

<p>Khain, Batygin &amp; Brown (2018) investigated the initial conditions required
to reproduce the observed distant Kuiper Belt. Starting from a narrow perihelion
distribution (q &#8712; 30–36 AU) produces only anti-aligned objects
(&#916;&#982; &#8776; &#960;), while a broad initial distribution (q &#8712; 30–300 AU)
generates both aligned and anti-aligned populations — matching observations. This
bimodality emerges naturally from the secular dynamics: anti-aligned particles are
confined by the quadrupole potential, while aligned particles occupy a secondary
stable island that exists only for moderate-to-high perihelion initial conditions.</p>

<h3>4.4 Mean-Motion Resonances</h3>

<p>Bailey, Brown &amp; Batygin (2018) examined whether the observed KBOs could be
in mean-motion resonance (MMR) with P9. Using Farey sequence enumeration to
catalog all possible N/M resonances, they found that the probability of the six
observed KBOs all residing in simple N/1 or N/2 resonances is less than 5%. While
individual resonances can constrain a<sub>9</sub>, the full resonance spectrum
produces a broad plateau rather than a sharp peak, limiting the diagnostic power
of resonance-based searches.</p>

<h3>4.5 Inner Oort Cloud Injection</h3>

<p>Batygin &amp; Brown (2021) showed that P9 can inject inner Oort Cloud (IOC)
objects — emplaced during the Sun's birth cluster phase (Plummer sphere with
M &#8776; 1200 M<sub>&#8857;</sub>, r &#8776; 0.35 pc) — into the distant Kuiper Belt.
IOC objects show weaker longitude-of-perihelion confinement (~67%) compared to
scattered disk objects (~88%), and preferentially populate the a &gt; 2000 AU
region, potentially contributing to the observed population of extreme TNOs.</p>
"#,
    );
}

fn write_section_5_statistical_evidence(html: &mut String) {
    html.push_str(
        r#"
<h2>5. Statistical Evidence</h2>

<h3>5.1 Accounting for Observational Bias (2017)</h3>

<p>A critical question is whether the observed clustering could be an artifact of
observational selection effects — surveys preferentially detect objects near
perihelion, where they are brightest, potentially creating spurious angular
concentrations. Brown (2017) addressed this with a detailed bias model
incorporating perihelion-distance detection probability (scaling as r<sup>&#8722;4</sup>
for flux-limited surveys) and galactic plane avoidance (reduced detection
efficiency at low galactic latitudes).</p>

<p>Using a sample of 10 KBOs with a &gt; 230 AU and a Monte Carlo framework with
10,000 bias-weighted random draws, the combined &#982; and &#969; clustering
significance was found to be 99.8% (p &#8776; 0.002) — comparable to the
unbiased result. The Rayleigh z-statistic, applied to the bias-corrected angular
distribution, confirmed that observational selection effects cannot account for
the observed clustering.</p>

<h3>5.2 Four-Dimensional Clustering (2019)</h3>

<p>Brown &amp; Batygin (2019) extended the analysis to 14 KBOs (adding
2013 RF<sub>98</sub>, 2014 SR<sub>349</sub>, 2015 GT<sub>50</sub>,
2015 KG<sub>163</sub>, 2015 RX<sub>245</sub>, and 2014 FE<sub>72</sub>)
and introduced Poincaré action-angle variables for a more rigorous treatment:</p>

<p>The canonical variables &#915; = 1 &#8722; &#8730;(1 &#8722; e<sup>2</sup>)
and Z = &#8730;(1 &#8722; e<sup>2</sup>)(1 &#8722; cos i) transform the clustering
analysis into a 4D problem: (x, y) = (&#8730;2&#915; cos &#982;, &#8730;2&#915; sin &#982;)
for perihelion direction, and (p, q) = (&#8730;2Z cos &#937;, &#8730;2Z sin &#937;)
for orbital pole position. The combined 4D clustering significance is 99.8%,
strengthening the case that the anomaly is physical rather than observational.</p>
"#,
    );
}

fn write_section_6_parameter_evolution(html: &mut String) {
    html.push_str(
        r#"
<h2>6. Parameter Evolution</h2>

<p>One of the most instructive aspects of the Planet Nine literature is the
progressive refinement of the predicted orbital parameters. The initial 2016
proposal was based on a single representative simulation; subsequent work
explored the full viable parameter space and incorporated additional constraints.</p>
"#,
    );

    // Embed the parameter evolution diagram
    let svg = diagrams::parameter_evolution_diagram();
    writeln!(
        html,
        r#"<div class="figure">{}<p class="figure-caption"><strong>Figure 1.</strong> Evolution of Planet Nine parameter estimates from 2016 to 2024. Points show best-fit or median values; vertical bars indicate 1&#963; ranges where available. The mass has decreased from 10 to ~6 M<sub>&#8853;</sub>, the semi-major axis from 700 to ~380 AU, the eccentricity from 0.6 to ~0.2, and the inclination from 30° to ~16°.</p></div>"#,
        svg
    )
    .unwrap();

    html.push_str(r#"
<h3>6.1 The 2019 Revision</h3>

<p>The most significant parameter update came from Batygin, Adams, Brown &amp; Becker
(2019), who conducted 1,794 simulations (1,134 semi-averaged + 660 fully resolved
N-body) spanning a systematic grid in (a<sub>9</sub>, e<sub>9</sub>, i<sub>9</sub>,
m<sub>9</sub>). The revised viable ranges — m<sub>9</sub> &#8776; 5–10 M<sub>&#8853;</sub>,
a<sub>9</sub> &#8776; 400–800 AU, e<sub>9</sub> &#8776; 0.2–0.5 — represented a
substantial shift from the 2016 values, driven primarily by the recognition that
lower-mass, lower-eccentricity orbits more robustly reproduce the observed KBO
population.</p>

<h3>6.2 The MCMC Posterior (2021)</h3>

<p>Brown &amp; Batygin (2021) constructed a formal statistical posterior using
Markov Chain Monte Carlo sampling, fitting to the orbital elements of 11 extreme
TNOs with a &gt; 250 AU. The resulting asymmetric Gaussian distributions yielded
the most precise parameter estimates to date:</p>

<table>
<tr><th>Parameter</th><th>2016</th><th>2019 Range</th><th>2021 MCMC</th></tr>
<tr><td>Mass (M<sub>&#8853;</sub>)</td><td>10</td><td>5–10</td><td>6.2<sup>+2.2</sup><sub>&#8722;1.3</sub></td></tr>
<tr><td>Semi-major axis (AU)</td><td>700</td><td>400–800</td><td>380<sup>+140</sup><sub>&#8722;80</sub></td></tr>
<tr><td>Eccentricity</td><td>0.6</td><td>0.2–0.5</td><td>~0.21</td></tr>
<tr><td>Inclination (°)</td><td>30</td><td>15–25</td><td>16 &pm; 5</td></tr>
<tr><td>Perihelion (AU)</td><td>280</td><td>200–400</td><td>300<sup>+85</sup><sub>&#8722;60</sub></td></tr>
<tr><td>V magnitude</td><td>~25</td><td>22–25</td><td>~22</td></tr>
</table>

<p>The convergence toward lower mass, smaller semi-major axis, and lower eccentricity
is physically motivated: a closer, less eccentric orbit produces stronger secular
coupling per unit mass, allowing a lighter planet to achieve the same dynamical
effect.</p>
"#);
}

fn write_section_7_observational_campaigns(html: &mut String) {
    html.push_str(
        r#"
<h2>7. Observational Search Campaigns</h2>

<p>Three major wide-field surveys have been systematically analyzed for Planet Nine
candidates, progressively shrinking the viable parameter space.</p>
"#,
    );

    // Embed the exclusion diagram
    let svg = diagrams::exclusion_timeline_diagram();
    writeln!(
        html,
        r#"<div class="figure">{}<p class="figure-caption"><strong>Figure 2.</strong> Cumulative exclusion of the Planet Nine parameter space by survey. ZTF (2021) excluded 56.4%, DES (2022) added 5.0% unique exclusion, and Pan-STARRS1 (2024) added 17.1%, for a combined total of 78.5%. The remaining 21.5% is concentrated at V &#8819; 22, requiring next-generation surveys.</p></div>"#,
        svg
    )
    .unwrap();

    html.push_str(
        r#"
<h3>7.1 Zwicky Transient Facility (2021)</h3>

<p>Brown &amp; Batygin (2021) searched the ZTF public archive (limiting magnitude
r &#8776; 20.5, sky coverage ~75%, tracklet linking efficiency 99.66%) against
100,000 synthetic P9 orbits drawn from the MCMC posterior. The detection pipeline
applied depth, footprint, and linking cuts with self-calibration corrections for
trailing losses and crowded fields. No candidates were identified, excluding
56.4% of the prior parameter space — primarily the brighter (V &lt; 20.5),
closer (a &lt; 400 AU) portion.</p>

<h3>7.2 Dark Energy Survey (2022)</h3>

<p>Belyakov, Bernardinelli &amp; Brown (2022) analyzed six years of DES data
(5,000 deg<sup>2</sup> footprint, 575 observing nights, 10&#963; depth g = 24.1).
Five surface color models were tested — fiducial, Neptune-like, methane at 40K,
super-Ganymede, and super-KBO — to account for uncertainty in P9's atmospheric
composition. Of synthetic orbits crossing the DES footprint, 87% were recovered.
DES contributes 5.0% unique exclusion (orbits not already excluded by ZTF),
raising the combined total to 61.2%.</p>

<h3>7.3 Pan-STARRS1 (2024)</h3>

<p>Brown, Holman &amp; Batygin (2024) completed the current survey trilogy with a
Pan-STARRS1 analysis (depth V &#8776; 21.5, linking threshold 9 detections,
efficiency 99.2%). PS1's complementary sky coverage and intermediate depth between
ZTF and DES proved particularly effective: 69,802 of 100,000 synthetic orbits were
detected, with 17,054 unique to PS1. This added 17.1% unique exclusion, bringing
the combined three-survey total to 78.5%.</p>

<p>The updated parameter estimates incorporating all three survey non-detections are:
a<sub>9</sub> = 500 &pm; 170 AU, m<sub>9</sub> = 6.6 &pm; 2.2 M<sub>&#8853;</sub>,
V = 22.0 &pm; 1.3 — shifted toward larger distance and fainter magnitude relative
to the unconstrained MCMC posterior.</p>
"#,
    );
}

fn write_section_8_latest_evidence(html: &mut String) {
    html.push_str(r#"
<h2>8. Latest Dynamical Evidence: Neptune-Crossing TNOs (2024)</h2>

<p>The most recent and arguably most compelling dynamical evidence for Planet Nine
comes from Batygin, Morbidelli, Brown &amp; Nesvorny (2024), who analyzed a sample
of 17 well-characterized Neptune-crossing TNOs satisfying a &gt; 100 AU,
i &lt; 40°, q &lt; 30 AU. These objects represent a distinct population from the
distant KBOs used in the original clustering argument — they interact directly
with Neptune and their orbital evolution is strongly influenced by close encounters.</p>

<p>The key insight is that in a P9-inclusive solar system, Neptune-crossing objects
develop a characteristic perihelion distance distribution shaped by P9's secular
perturbations. The authors compared two models:</p>

<table>
<tr><th>Model</th><th>&#950; Statistic</th><th>p-value</th><th>Significance</th></tr>
<tr><td>P9-inclusive (m=5, a=500, e=0.25, i=20°)</td><td>&#8722;7.9</td><td>0.41</td><td>Consistent</td></tr>
<tr><td>P9-free null hypothesis</td><td>&#8722;16.5</td><td>0.0034</td><td>~5&#963; rejection</td></tr>
</table>

<p>The P9-free model is rejected at approximately 5&#963; significance via the
Kolmogorov-Smirnov test on perihelion distributions, while the P9-inclusive model
is fully consistent with observations (p = 0.41). This result is particularly
powerful because it uses an independent population (Neptune-crossers rather than
distant KBOs), an independent observable (perihelion distribution rather than
angular clustering), and cannot be attributed to observational bias (Neptune-crossers
are detected at perihelion regardless of angular position).</p>
"#);
}

fn write_section_9_current_state(html: &mut String) {
    html.push_str(
        r#"
<h2>9. Current State of the Hypothesis</h2>
"#,
    );

    // Embed the synthesis diagram
    let svg = diagrams::hypothesis_synthesis_diagram();
    writeln!(
        html,
        r#"<div class="figure">{}<p class="figure-caption"><strong>Figure 3.</strong> Synthesis of the Planet Nine hypothesis across all papers. Each row represents a paper's contribution, organized by type: dynamical evidence (blue), statistical analysis (green), parameter refinement (gold), and observational searches (red). The hypothesis has accumulated multiple independent lines of evidence while the predicted parameter space has been progressively narrowed.</p></div>"#,
        svg
    )
    .unwrap();

    html.push_str(
        r#"
<h3>9.1 Lines of Evidence</h3>

<p>The Planet Nine hypothesis is supported by multiple independent lines of evidence,
each addressing a distinct dynamical phenomenon:</p>

<div class="highlight">
<strong>Six independent signatures explained by a single body:</strong>
<ol style="margin-top:0.5rem; padding-left:1.5rem;">
<li>Longitude of perihelion clustering of distant KBOs (&#982; &#8776; 71°)</li>
<li>Solar obliquity (6° tilt from invariable plane)</li>
<li>Highly inclined and retrograde TNOs (Drac, Niku, 2016 NM<sub>56</sub>)</li>
<li>Bimodal distant Kuiper Belt (aligned + anti-aligned populations)</li>
<li>Inner Oort Cloud injection into the distant Kuiper Belt</li>
<li>Neptune-crossing TNO perihelion distribution (~5&#963; rejection of P9-free)</li>
</ol>
</div>

<h3>9.2 Remaining Parameter Space</h3>

<p>As of 2024, the combined ZTF + DES + Pan-STARRS1 non-detections have excluded
78.5% of the viable parameter space. The remaining 21.5% is characterized by:</p>

<table>
<tr><th>Parameter</th><th>Remaining Range</th><th>Peak Probability</th></tr>
<tr><td>Semi-major axis</td><td>400–700 AU</td><td>~500 AU</td></tr>
<tr><td>Mass</td><td>4–9 M<sub>&#8853;</sub></td><td>~6.5 M<sub>&#8853;</sub></td></tr>
<tr><td>V magnitude</td><td>21.5–24</td><td>~22</td></tr>
<tr><td>Eccentricity</td><td>0.15–0.35</td><td>~0.2</td></tr>
<tr><td>Inclination</td><td>10–25°</td><td>~16°</td></tr>
</table>

<p>The remaining orbits are concentrated at fainter magnitudes (V &#8819; 22) where
current surveys lack sufficient depth or coverage. The object is most likely near
aphelion in a region of sky not yet surveyed to adequate depth.</p>
"#,
    );
}

fn write_section_10_future(html: &mut String) {
    html.push_str(
        r#"
<h2>10. Future Prospects</h2>

<h3>10.1 Vera C. Rubin Observatory (LSST)</h3>

<p>The Legacy Survey of Space and Time (LSST) at the Vera C. Rubin Observatory
represents the most promising near-term opportunity for Planet Nine detection.
With a single-visit depth of r &#8776; 24.5 and co-added depth exceeding 27th
magnitude over 18,000 deg<sup>2</sup>, LSST will probe essentially the entire
remaining parameter space. At the current best-estimate distance of ~500 AU and
magnitude V &#8776; 22, P9 would be well within LSST's detection capability.</p>

<p>First light is expected in 2025, with the full survey commencing shortly
thereafter. Within the first 1–2 years of operation, LSST should have sufficient
temporal baseline for motion-based identification of a slow-moving outer solar
system object at P9's expected angular velocity (~1–3 arcsec/hour).</p>

<h3>10.2 Theoretical Refinements</h3>

<p>Several theoretical directions could further constrain or support the hypothesis:</p>

<p><strong>Discovery of new extreme TNOs.</strong> Each additional distant KBO with
a &gt; 250 AU provides additional clustering data. The current sample of ~14–17
objects (depending on selection criteria) produces 99.8% confidence in clustering;
doubling this sample would either dramatically strengthen or weaken the case.</p>

<p><strong>Improved bias modeling.</strong> More sophisticated treatment of survey
selection functions, particularly for the DES and forthcoming LSST data, could
sharpen the statistical significance of clustering measurements.</p>

<p><strong>Alternative explanations.</strong> The "Planet Nine-free" hypothesis
continues to face the challenge of simultaneously explaining all six dynamical
signatures listed in §9.1. Proposed alternatives — a massive primordial disk,
a rogue planet flyby, or self-gravity of the distant Kuiper Belt — each address
at most one or two of these phenomena.</p>

<h3>10.3 What Detection Would Tell Us</h3>

<p>If confirmed, Planet Nine would be the first new planet discovered in the solar
system since Neptune in 1846. At 5–7 M<sub>&#8853;</sub>, it would be a
"super-Earth" or "sub-Neptune" — the most common type of planet in the galaxy
according to exoplanet surveys, yet conspicuously absent from our own solar
system. Its formation and emplacement mechanism (in situ formation, scattering
from the giant planet region, or capture from a passing star) would carry
profound implications for solar system formation models.</p>

<p>The convergence of multiple independent dynamical signatures, strengthening
statistical evidence (now at ~5&#963; from Neptune-crossing TNOs alone), and
the systematic exclusion of 78% of parameter space paints a picture of a
hypothesis approaching its moment of truth. The remaining 22% of parameter
space is squarely within reach of next-generation surveys. Within the decade,
we will likely know whether the solar system's most distant planet is real.</p>
"#,
    );
}

fn write_references(html: &mut String) {
    html.push_str(r#"
<h2>References</h2>

<div class="ref-list">
<p>[1] Batygin, K. &amp; Brown, M. E. (2016a). "Evidence for a Distant Giant Planet in the Solar System." <em>The Astronomical Journal</em>, 151(2), 22.</p>
<p>[2] Brown, M. E. &amp; Batygin, K. (2016). "Observational Constraints on the Orbit and Location of Planet Nine in the Outer Solar System." <em>The Astrophysical Journal Letters</em>, 824(2), L23.</p>
<p>[3] Bailey, E., Batygin, K. &amp; Brown, M. E. (2016). "Solar Obliquity Induced by Planet Nine." <em>The Astronomical Journal</em>, 152(5), 126.</p>
<p>[4] Batygin, K. &amp; Brown, M. E. (2016b). "Generation of Highly Inclined Trans-Neptunian Objects by Planet Nine." <em>The Astrophysical Journal Letters</em>, 833(1), L3.</p>
<p>[5] Brown, M. E. (2017). "Observational Bias and the Clustering of Distant Eccentric Kuiper Belt Objects." <em>The Astronomical Journal</em>, 154(2), 65.</p>
<p>[6] Khain, T., Batygin, K. &amp; Brown, M. E. (2018). "The Generation of the Distant Kuiper Belt by Planet Nine from an Initially Broad Perihelion Distribution." <em>The Astronomical Journal</em>, 155(6), 250.</p>
<p>[7] Bailey, E., Brown, M. E. &amp; Batygin, K. (2018). "Feasibility of a Resonance-Based Planet Nine Search." <em>The Astronomical Journal</em>, 156(2), 74.</p>
<p>[8] Brown, M. E. &amp; Batygin, K. (2019). "Orbital Clustering in the Distant Solar System." <em>The Astronomical Journal</em>, 157(2), 62.</p>
<p>[9] Batygin, K., Adams, F. C., Brown, M. E. &amp; Becker, J. C. (2019). "The Planet Nine Hypothesis." <em>Physics Reports</em>, 805, 1–53.</p>
<p>[10] Batygin, K. &amp; Brown, M. E. (2021). "Injection of Inner Oort Cloud Objects into the Distant Kuiper Belt by Planet Nine." <em>The Astrophysical Journal Letters</em>, 910(2), L20.</p>
<p>[11] Brown, M. E. &amp; Batygin, K. (2021). "The Orbit of Planet Nine." <em>The Astronomical Journal</em>, 162(5), 219.</p>
<p>[12] Brown, M. E. &amp; Batygin, K. (2021). "A Search for Planet Nine Using the Zwicky Transient Facility." <em>The Astronomical Journal</em>, 163(3), 102.</p>
<p>[13] Belyakov, M., Bernardinelli, P. H. &amp; Brown, M. E. (2022). "Limits on the Detection of Planet Nine in the Dark Energy Survey." <em>The Planetary Science Journal</em>, 3(7), 159.</p>
<p>[14] Brown, M. E., Holman, M. J. &amp; Batygin, K. (2024). "A Pan-STARRS1 Search for Planet Nine." <em>The Astronomical Journal</em>.</p>
<p>[15] Batygin, K., Morbidelli, A., Brown, M. E. &amp; Nesvorny, D. (2024). "Generation of Low-Inclination, Neptune-Crossing Trans-Neptunian Objects by Planet Nine." <em>The Astrophysical Journal Letters</em>.</p>
</div>
"#);
}

fn write_footer(html: &mut String) {
    html.push_str(
        r#"
</body>
</html>
"#,
    );
}
