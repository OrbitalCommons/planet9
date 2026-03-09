use std::fmt::Write;

// Helper to avoid raw string delimiter conflicts with SVG color codes (#xxx).
// All SVG write macros use r##"..."## delimiters.

/// Figure 1: Parameter evolution from 2016 to 2024
/// Shows four subplots: mass, semi-major axis, eccentricity, inclination
pub fn parameter_evolution_diagram() -> String {
    let mut svg = String::with_capacity(16_000);
    let w = 800;
    let h = 600;

    writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Georgia, serif">"##,
        w, h
    )
    .unwrap();

    writeln!(
        svg,
        r##"<rect width="{}" height="{}" fill="white" rx="4"/>"##,
        w, h
    )
    .unwrap();

    writeln!(
        svg,
        r##"<text x="{}" y="25" text-anchor="middle" font-size="14" font-weight="bold" fill="#2c5f8a">Evolution of Planet Nine Parameter Estimates (2016–2024)</text>"##,
        w / 2
    )
    .unwrap();

    // Four subplots in a 2x2 grid
    // (ox, oy, title, unit, data: [(year_label, year_x, val, err)], y_min, y_max)
    let plots: [(f64, f64, &str, &str, &[(&str, f64, f64, f64)], f64, f64); 4] = [
        (
            60.0,
            50.0,
            "Mass",
            "M\u{2295}",
            &[
                ("2016", 0.0, 10.0, 0.0),
                ("2019", 0.4, 7.5, 2.5),
                ("2021", 0.7, 6.2, 2.2),
                ("2024", 1.0, 6.6, 2.2),
            ],
            0.0,
            15.0,
        ),
        (
            440.0,
            50.0,
            "Semi-major Axis",
            "AU",
            &[
                ("2016", 0.0, 700.0, 0.0),
                ("2019", 0.4, 600.0, 200.0),
                ("2021", 0.7, 380.0, 140.0),
                ("2024", 1.0, 500.0, 170.0),
            ],
            100.0,
            1000.0,
        ),
        (
            60.0,
            330.0,
            "Eccentricity",
            "",
            &[
                ("2016", 0.0, 0.6, 0.0),
                ("2019", 0.4, 0.35, 0.15),
                ("2021", 0.7, 0.21, 0.10),
                ("2024", 1.0, 0.20, 0.10),
            ],
            0.0,
            0.9,
        ),
        (
            440.0,
            330.0,
            "Inclination",
            "\u{00B0}",
            &[
                ("2016", 0.0, 30.0, 0.0),
                ("2019", 0.4, 20.0, 5.0),
                ("2021", 0.7, 16.0, 5.0),
                ("2024", 1.0, 16.0, 5.0),
            ],
            0.0,
            40.0,
        ),
    ];

    let plot_w = 300.0;
    let plot_h = 220.0;
    let colors = ["#d32f2f", "#e65100", "#1565c0", "#2e7d32"];

    for (ox, oy, title, unit, data, y_min, y_max) in &plots {
        let px = ox + 50.0;
        let py = oy + 30.0;
        let pw = plot_w - 70.0;
        let ph = plot_h - 60.0;
        let y_range = y_max - y_min;

        // Axes background
        writeln!(
            svg,
            r##"<rect x="{px}" y="{py}" width="{pw}" height="{ph}" fill="#f8f8f8" stroke="#ccc"/>"##
        )
        .unwrap();

        // Title
        writeln!(svg, r##"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-weight="bold" fill="#333">{title}</text>"##,
            px + pw / 2.0, *oy + 22.0).unwrap();

        // Y-axis ticks and grid
        let n_ticks = 5;
        for i in 0..=n_ticks {
            let frac = i as f64 / n_ticks as f64;
            let val = y_min + frac * y_range;
            let y_pos = py + ph - frac * ph;

            writeln!(svg, r##"<line x1="{px}" y1="{y_pos}" x2="{}" y2="{y_pos}" stroke="#eee" stroke-width="0.5"/>"##, px + pw).unwrap();

            let label = if y_range < 1.0 {
                format!("{:.1}", val)
            } else {
                format!("{:.0}", val)
            };
            writeln!(svg, r##"<text x="{}" y="{}" text-anchor="end" font-size="9" fill="#666">{label}</text>"##,
                px - 4.0, y_pos + 3.0).unwrap();
        }

        // Unit label
        if !unit.is_empty() {
            writeln!(svg, r##"<text x="{}" y="{}" text-anchor="end" font-size="9" fill="#999">({unit})</text>"##,
                px - 4.0, py - 4.0).unwrap();
        }

        // Data points with error bars and connecting line
        let mut path = String::new();
        for (idx, (year_label, year_x, val, err)) in data.iter().enumerate() {
            let x = px + 20.0 + year_x * (pw - 40.0);
            let y = py + ph - (val - y_min) / y_range * ph;
            let color = colors[idx];

            // Error bar
            if *err > 0.0 {
                let y_top = (py + ph - ((val + err) - y_min) / y_range * ph).max(py);
                let y_bot = (py + ph - ((val - err) - y_min) / y_range * ph).min(py + ph);
                writeln!(svg, r##"<line x1="{x}" y1="{y_top}" x2="{x}" y2="{y_bot}" stroke="{color}" stroke-width="2" opacity="0.4"/>"##).unwrap();
                writeln!(svg, r##"<line x1="{}" y1="{y_top}" x2="{}" y2="{y_top}" stroke="{color}" stroke-width="2" opacity="0.4"/>"##, x - 4.0, x + 4.0).unwrap();
                writeln!(svg, r##"<line x1="{}" y1="{y_bot}" x2="{}" y2="{y_bot}" stroke="{color}" stroke-width="2" opacity="0.4"/>"##, x - 4.0, x + 4.0).unwrap();
            }

            // Point
            writeln!(svg, r##"<circle cx="{x}" cy="{y}" r="5" fill="{color}" stroke="white" stroke-width="1.5"/>"##).unwrap();

            // Year label
            writeln!(svg, r##"<text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="#333">{year_label}</text>"##, py + ph + 14.0).unwrap();

            // Connecting line path
            if idx == 0 {
                write!(path, "M{x},{y}").unwrap();
            } else {
                write!(path, " L{x},{y}").unwrap();
            }
        }

        writeln!(svg, r##"<path d="{path}" fill="none" stroke="#2c5f8a" stroke-width="1.5" stroke-dasharray="4,3" opacity="0.5"/>"##).unwrap();
    }

    svg.push_str("</svg>");
    svg
}

/// Figure 2: Survey exclusion timeline showing cumulative parameter space reduction
pub fn exclusion_timeline_diagram() -> String {
    let mut svg = String::with_capacity(12_000);
    let w = 800;
    let h = 400;

    writeln!(svg, r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Georgia, serif">"##).unwrap();
    writeln!(
        svg,
        r##"<rect width="{w}" height="{h}" fill="white" rx="4"/>"##
    )
    .unwrap();
    writeln!(svg, r##"<text x="{}" y="28" text-anchor="middle" font-size="14" font-weight="bold" fill="#2c5f8a">Cumulative Exclusion of Planet Nine Parameter Space</text>"##, w / 2).unwrap();

    let ox = 120.0_f64;
    let oy = 50.0_f64;
    let pw = 560.0_f64;
    let ph = 280.0_f64;

    writeln!(
        svg,
        r##"<rect x="{ox}" y="{oy}" width="{pw}" height="{ph}" fill="#f8f8f8" stroke="#ccc"/>"##
    )
    .unwrap();

    // Y-axis (percentage)
    for i in 0..=10 {
        let pct = i as f64 * 10.0;
        let y = oy + ph - pct / 100.0 * ph;
        writeln!(
            svg,
            r##"<line x1="{ox}" y1="{y}" x2="{}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>"##,
            ox + pw
        )
        .unwrap();
        writeln!(
            svg,
            r##"<text x="{}" y="{}" text-anchor="end" font-size="10" fill="#666">{:.0}%</text>"##,
            ox - 6.0,
            y + 3.5,
            pct
        )
        .unwrap();
    }

    // Y-axis title
    writeln!(svg, r##"<text x="15" y="{}" text-anchor="middle" font-size="11" fill="#333" transform="rotate(-90, 15, {})">Parameter Space Excluded</text>"##,
        oy + ph / 2.0, oy + ph / 2.0).unwrap();

    // Milestones: (x_frac, unique_pct, cumulative_pct, label, year, color)
    let milestones: [(f64, f64, f64, &str, &str, &str); 5] = [
        (0.0, 0.0, 0.0, "Pre-survey", "2016", "#999"),
        (0.3, 56.4, 56.4, "ZTF", "2021", "#e53935"),
        (0.55, 5.0, 61.2, "DES", "2022", "#fb8c00"),
        (0.8, 17.1, 78.5, "PS1", "2024", "#1e88e5"),
        (1.0, 0.0, 78.5, "Current", "2024+", "#43a047"),
    ];

    // Filled area under cumulative curve
    let mut area_path = format!("M{ox},{}", oy + ph);
    for &(x_frac, _, cum_pct, _, _, _) in &milestones {
        let x = ox + x_frac * pw;
        let y = oy + ph - cum_pct / 100.0 * ph;
        write!(area_path, " L{x},{y}").unwrap();
    }
    write!(area_path, " L{},{} Z", ox + pw, oy + ph).unwrap();
    writeln!(
        svg,
        r##"<path d="{area_path}" fill="#2c5f8a" opacity="0.12"/>"##
    )
    .unwrap();

    // Remaining space area (top)
    let mut remaining_path = format!("M{ox},{oy}");
    for &(x_frac, _, cum_pct, _, _, _) in &milestones {
        let x = ox + x_frac * pw;
        let y = oy + ph - cum_pct / 100.0 * ph;
        write!(remaining_path, " L{x},{y}").unwrap();
    }
    write!(remaining_path, " L{},{oy} Z", ox + pw).unwrap();
    writeln!(
        svg,
        r##"<path d="{remaining_path}" fill="#43a047" opacity="0.08"/>"##
    )
    .unwrap();

    // Cumulative line
    let mut line_path = String::new();
    for (i, &(x_frac, _, cum_pct, _, _, _)) in milestones.iter().enumerate() {
        let x = ox + x_frac * pw;
        let y = oy + ph - cum_pct / 100.0 * ph;
        if i == 0 {
            write!(line_path, "M{x},{y}").unwrap();
        } else {
            write!(line_path, " L{x},{y}").unwrap();
        }
    }
    writeln!(
        svg,
        r##"<path d="{line_path}" fill="none" stroke="#2c5f8a" stroke-width="2.5"/>"##
    )
    .unwrap();

    // Points and labels
    for &(x_frac, unique_pct, cum_pct, label, year, color) in &milestones {
        let x = ox + x_frac * pw;
        let y = oy + ph - cum_pct / 100.0 * ph;

        writeln!(
            svg,
            r##"<circle cx="{x}" cy="{y}" r="6" fill="{color}" stroke="white" stroke-width="2"/>"##
        )
        .unwrap();

        let label_y = if cum_pct > 50.0 { y - 16.0 } else { y - 14.0 };
        writeln!(svg, r##"<text x="{x}" y="{label_y}" text-anchor="middle" font-size="11" font-weight="bold" fill="{color}">{label}</text>"##).unwrap();
        writeln!(svg, r##"<text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="#666">{year}</text>"##, oy + ph + 16.0).unwrap();

        if unique_pct > 0.0 {
            writeln!(svg, r##"<text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="{color}">+{unique_pct:.1}%</text>"##, label_y - 13.0).unwrap();
        }
    }

    // Region labels
    writeln!(svg, r##"<text x="{}" y="{}" text-anchor="middle" font-size="13" fill="#43a047" font-weight="bold">Remaining: 21.5%</text>"##,
        ox + pw * 0.75, oy + 30.0).unwrap();
    writeln!(svg, r##"<text x="{}" y="{}" text-anchor="middle" font-size="13" fill="#2c5f8a" font-weight="bold">Excluded: 78.5%</text>"##,
        ox + pw * 0.75, oy + ph - 20.0).unwrap();

    // LSST projected dashed extension
    let last_x = ox + pw;
    let last_y = oy + ph - 78.5 / 100.0 * ph;
    let lsst_y = oy + ph - 95.0 / 100.0 * ph;
    writeln!(svg, r##"<line x1="{}" y1="{last_y}" x2="{}" y2="{lsst_y}" stroke="#43a047" stroke-width="1.5" stroke-dasharray="6,4"/>"##,
        last_x - 10.0, last_x - 20.0).unwrap();
    writeln!(svg, r##"<text x="{}" y="{}" text-anchor="start" font-size="9" fill="#43a047" font-style="italic">LSST?</text>"##,
        ox + pw - 50.0, lsst_y + 4.0).unwrap();

    svg.push_str("</svg>");
    svg
}

/// Figure 3: Comprehensive hypothesis synthesis — a timeline/contribution matrix
pub fn hypothesis_synthesis_diagram() -> String {
    let mut svg = String::with_capacity(24_000);
    let w = 850;
    let h = 720;

    writeln!(svg, r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Georgia, serif">"##).unwrap();
    writeln!(
        svg,
        r##"<rect width="{w}" height="{h}" fill="white" rx="4"/>"##
    )
    .unwrap();

    // Title
    writeln!(svg, r##"<text x="{}" y="25" text-anchor="middle" font-size="14" font-weight="bold" fill="#2c5f8a">Planet Nine Hypothesis: Paper Contributions &amp; Evidence Map</text>"##, w / 2).unwrap();

    // Arrow marker definition (put defs early)
    writeln!(svg, r##"<defs><marker id="arrowhead" markerWidth="10" markerHeight="7" refX="10" refY="3.5" orient="auto"><polygon points="0 0, 10 3.5, 0 7" fill="#2c5f8a"/></marker></defs>"##).unwrap();

    // Legend
    let legend_items: [(&str, &str); 5] = [
        ("#1565c0", "Dynamical Evidence"),
        ("#2e7d32", "Statistical Analysis"),
        ("#e65100", "Parameter Refinement"),
        ("#c62828", "Observational Search"),
        ("#6a1b9a", "Review / Synthesis"),
    ];
    let legend_y = 42.0_f64;
    for (i, (color, label)) in legend_items.iter().enumerate() {
        let x = 100.0 + i as f64 * 150.0;
        writeln!(
            svg,
            r##"<rect x="{x}" y="{legend_y}" width="12" height="12" fill="{color}" rx="2"/>"##
        )
        .unwrap();
        writeln!(
            svg,
            r##"<text x="{}" y="{}" font-size="9" fill="#333">{label}</text>"##,
            x + 16.0,
            legend_y + 10.0
        )
        .unwrap();
    }

    // Paper data
    let papers: Vec<(&str, &str, &str, &str, [&str; 3], &str)> = vec![
        (
            "2016",
            "Evidence for P9",
            "Batygin &amp; Brown",
            "#1565c0",
            [
                "6 KBO clustering",
                "Secular confinement",
                "Anti-aligned orbits",
            ],
            "p = 0.007%",
        ),
        (
            "2016",
            "Constraints",
            "Brown &amp; Batygin",
            "#e65100",
            [
                "19\u{00D7}9\u{00D7}5 param grid",
                "V ~ 22\u{2013}25 predicted",
                "Perihelion band",
            ],
            "a = 200\u{2013}2000",
        ),
        (
            "2016",
            "Solar Obliquity",
            "Bailey+",
            "#1565c0",
            [
                "6\u{00B0} tilt explained",
                "Spin-orbit coupling",
                "i\u{2089} ~ 15\u{2013}25\u{00B0}",
            ],
            "6\u{00B0} obliquity",
        ),
        (
            "2016",
            "Inclined TNOs",
            "Batygin &amp; Brown",
            "#1565c0",
            [
                "Kozai-Lidov pump",
                "Retrograde orbits",
                "Drac/Niku explained",
            ],
            "i &gt; 90\u{00B0}",
        ),
        (
            "2017",
            "Bias Analysis",
            "Brown",
            "#2e7d32",
            [
                "r\u{207B}\u{2074} bias model",
                "MC rejection test",
                "Bias-corrected sig.",
            ],
            "99.8% conf.",
        ),
        (
            "2018",
            "Kuiper Belt Gen.",
            "Khain+",
            "#1565c0",
            [
                "Broad q\u{2080} required",
                "Bimodal population",
                "Aligned + anti-aligned",
            ],
            "2 populations",
        ),
        (
            "2018",
            "Resonance Search",
            "Bailey+",
            "#e65100",
            [
                "Farey seq. catalog",
                "N/1, N/2 resonances",
                "Broad a\u{2089} plateau",
            ],
            "P(res) &lt; 5%",
        ),
        (
            "2019",
            "4D Clustering",
            "Brown &amp; Batygin",
            "#2e7d32",
            [
                "14 KBO sample",
                "Poincar\u{00E9} variables",
                "\u{03D6} + pole combined",
            ],
            "99.8% (4D)",
        ),
        (
            "2019",
            "Review",
            "Batygin+",
            "#6a1b9a",
            [
                "1794 simulations",
                "Params revised",
                "m = 5\u{2013}10 M\u{2295}",
            ],
            "a = 400\u{2013}800",
        ),
        (
            "2021",
            "Oort Cloud Inj.",
            "Batygin &amp; Brown",
            "#1565c0",
            [
                "IOC \u{2192} KB injection",
                "67% vs 88% confine",
                "a &gt; 2000 AU pop.",
            ],
            "f = 67%",
        ),
        (
            "2021",
            "Orbit MCMC",
            "Brown &amp; Batygin",
            "#e65100",
            ["MCMC posterior", "11 ETNO fit", "Asymm. Gaussian"],
            "m = 6.2 M\u{2295}",
        ),
        (
            "2021",
            "ZTF Search",
            "Brown &amp; Batygin",
            "#c62828",
            ["V &lt; 20.5 searched", "99.66% link eff.", "56.4% excluded"],
            "\u{2013}56.4%",
        ),
        (
            "2022",
            "DES Limits",
            "Belyakov+",
            "#c62828",
            ["5 color models", "87% recovery", "+5% unique excl."],
            "\u{2013}61.2%",
        ),
        (
            "2024",
            "Pan-STARRS1",
            "Brown+",
            "#c62828",
            ["V &lt; 21.5 searched", "+17.1% unique", "78% total excl."],
            "\u{2013}78.5%",
        ),
        (
            "2024",
            "Neptune-crossing",
            "Batygin+",
            "#1565c0",
            ["17 N-crossing TNOs", "KS + A-D tests", "P9-free rejected"],
            "~5\u{03C3} rej.",
        ),
    ];

    let start_y = 70.0_f64;
    let row_h = 40.0_f64;
    let year_x = 30.0_f64;
    let label_x = 80.0_f64;
    let ref_x = 210.0_f64;
    let findings_x = 340.0_f64;
    let number_x = 740.0_f64;

    // Column headers
    for (x, text) in [
        (year_x, "Year"),
        (label_x, "Paper"),
        (ref_x, "Authors"),
        (findings_x + 80.0, "Key Findings"),
        (number_x, "Key Result"),
    ] {
        writeln!(
            svg,
            r##"<text x="{x}" y="{}" font-size="10" font-weight="bold" fill="#333">{text}</text>"##,
            start_y + 4.0
        )
        .unwrap();
    }

    // Header underline
    writeln!(
        svg,
        r##"<line x1="20" y1="{}" x2="{}" y2="{}" stroke="#2c5f8a" stroke-width="2"/>"##,
        start_y + 10.0,
        w - 20,
        start_y + 10.0
    )
    .unwrap();

    for (i, (year, paper_label, authors, cat_color, findings, key_num)) in papers.iter().enumerate()
    {
        let y = start_y + 14.0 + (i + 1) as f64 * row_h;

        // Alternating row background
        if i % 2 == 0 {
            writeln!(
                svg,
                r##"<rect x="20" y="{}" width="{}" height="{}" fill="#f5f5f5" rx="2"/>"##,
                y - 12.0,
                w - 40,
                row_h - 2.0
            )
            .unwrap();
        }

        // Category color bar
        writeln!(
            svg,
            r##"<rect x="20" y="{}" width="4" height="{}" fill="{cat_color}" rx="1"/>"##,
            y - 12.0,
            row_h - 2.0
        )
        .unwrap();

        // Year
        writeln!(
            svg,
            r##"<text x="{year_x}" y="{}" font-size="10" fill="#666">{year}</text>"##,
            y + 4.0
        )
        .unwrap();

        // Paper label
        writeln!(svg, r##"<text x="{label_x}" y="{}" font-size="10" font-weight="bold" fill="{cat_color}">{paper_label}</text>"##, y + 4.0).unwrap();

        // Authors
        writeln!(svg, r##"<text x="{ref_x}" y="{}" font-size="9" fill="#666" font-style="italic">{authors}</text>"##, y + 4.0).unwrap();

        // Findings
        for (j, finding) in findings.iter().enumerate() {
            let fx = findings_x + j as f64 * 130.0;
            writeln!(
                svg,
                r##"<text x="{fx}" y="{}" font-size="8.5" fill="#444">{}{finding}</text>"##,
                y + 4.0,
                "\u{2022} "
            )
            .unwrap();
        }

        // Key number
        writeln!(svg, r##"<text x="{number_x}" y="{}" font-size="10" font-weight="bold" fill="{cat_color}">{key_num}</text>"##, y + 4.0).unwrap();
    }

    // Timeline arrow at bottom
    let arrow_y = start_y + 14.0 + (papers.len() + 1) as f64 * row_h;
    writeln!(svg, r##"<line x1="30" y1="{arrow_y}" x2="{}" y2="{arrow_y}" stroke="#2c5f8a" stroke-width="2" marker-end="url(#arrowhead)"/>"##, w - 40).unwrap();

    // Timeline phase labels
    let timeline_labels: [(f64, &str, &str); 5] = [
        (0.05, "Initial", "Proposal"),
        (0.25, "Dynamical", "Consequences"),
        (0.50, "Statistical", "Refinement"),
        (0.75, "Survey", "Exclusion"),
        (0.95, "Strongest", "Evidence"),
    ];
    for (frac, line1, line2) in &timeline_labels {
        let x = 30.0 + frac * (w as f64 - 70.0);
        writeln!(svg, r##"<text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="#2c5f8a" font-weight="bold">{line1}</text>"##, arrow_y + 14.0).unwrap();
        writeln!(svg, r##"<text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="#2c5f8a">{line2}</text>"##, arrow_y + 25.0).unwrap();
    }

    svg.push_str("</svg>");
    svg
}
