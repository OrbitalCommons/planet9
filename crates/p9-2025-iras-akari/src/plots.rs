//! SVG plot generation for the IRAS-AKARI Planet Nine search.
//!
//! Produces:
//! - Flux detectability plot: predicted flux vs distance for IRAS/AKARI limits
//! - Proper motion diagram: angular separation vs distance for the 23-yr baseline
//! - Search window visualization in (separation, flux ratio) space

use std::fmt::Write;

/// Generate an SVG plot showing predicted thermal flux vs. heliocentric
/// distance, with IRAS and AKARI sensitivity limits overlaid.
///
/// `distances`: heliocentric distances in AU
/// `flux_60um`: predicted flux at 60 µm in Jy (one per distance)
/// `flux_90um`: predicted flux at 90 µm in Jy (one per distance)
pub fn flux_vs_distance_plot(
    distances: &[f64],
    flux_60um: &[f64],
    flux_90um: &[f64],
    width: u32,
    height: u32,
) -> String {
    assert_eq!(distances.len(), flux_60um.len());
    assert_eq!(distances.len(), flux_90um.len());
    if distances.is_empty() {
        return String::new();
    }

    let margin = 70.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let d_min = distances.first().copied().unwrap_or(200.0);
    let d_max = distances.last().copied().unwrap_or(1000.0);
    let d_range = d_max - d_min;

    // Use log scale for flux
    let all_flux: Vec<f64> = flux_60um.iter().chain(flux_90um.iter()).copied().collect();
    let f_min_log = all_flux
        .iter()
        .filter(|f| **f > 0.0)
        .map(|f| f.log10())
        .fold(f64::INFINITY, f64::min)
        - 0.5;
    let f_max_log = all_flux
        .iter()
        .filter(|f| **f > 0.0)
        .map(|f| f.log10())
        .fold(f64::NEG_INFINITY, f64::max)
        + 0.5;
    let f_range_log = f_max_log - f_min_log;

    let mut svg = String::with_capacity(4000);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">P9 Thermal Flux vs Distance</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();
    writeln!(
        svg,
        r##"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="#f8f8f8" stroke="black" stroke-width="1"/>"##,
    )
    .unwrap();

    // IRAS sensitivity line (0.2 Jy)
    let iras_y = margin + plot_h - ((0.2_f64.log10() - f_min_log) / f_range_log * plot_h);
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{iras_y:.1}" x2="{}" y2="{iras_y:.1}" stroke="orangered" stroke-width="1.5" stroke-dasharray="6,3"/>"#,
        margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif" fill="orangered">IRAS 60µm limit (0.2 Jy)</text>"#,
        margin + 5.0,
        iras_y - 4.0,
    )
    .unwrap();

    // AKARI sensitivity line (0.55 Jy)
    let akari_y = margin + plot_h - ((0.55_f64.log10() - f_min_log) / f_range_log * plot_h);
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{akari_y:.1}" x2="{}" y2="{akari_y:.1}" stroke="steelblue" stroke-width="1.5" stroke-dasharray="6,3"/>"#,
        margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif" fill="steelblue">AKARI 90µm limit (0.55 Jy)</text>"#,
        margin + 5.0,
        akari_y - 4.0,
    )
    .unwrap();

    // 60 µm curve
    draw_curve(
        &mut svg,
        distances,
        flux_60um,
        d_min,
        d_range,
        f_min_log,
        f_range_log,
        margin,
        plot_w,
        plot_h,
        "orangered",
    );
    // 90 µm curve
    draw_curve(
        &mut svg,
        distances,
        flux_90um,
        d_min,
        d_range,
        f_min_log,
        f_range_log,
        margin,
        plot_w,
        plot_h,
        "steelblue",
    );

    // Search window shading (500–700 AU)
    let x_500 = margin + (500.0 - d_min) / d_range * plot_w;
    let x_700 = margin + (700.0 - d_min) / d_range * plot_w;
    writeln!(
        svg,
        r#"<rect x="{x_500:.1}" y="{margin}" width="{:.1}" height="{plot_h}" fill="gold" fill-opacity="0.15" stroke="goldenrod" stroke-width="1" stroke-dasharray="3,3"/>"#,
        x_700 - x_500,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{}" font-size="9" font-family="sans-serif" fill="goldenrod" text-anchor="middle">Search Window</text>"#,
        (x_500 + x_700) / 2.0,
        margin + 15.0,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Heliocentric Distance (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 40.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Flux Density (Jy)</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // X-axis ticks
    for &tick in &[200.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0] {
        if tick >= d_min && tick <= d_max {
            let x = margin + (tick - d_min) / d_range * plot_w;
            writeln!(
                svg,
                r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="8" font-family="sans-serif">{tick:.0}</text>"#,
                margin + plot_h + 15.0,
            )
            .unwrap();
        }
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

fn draw_curve(
    svg: &mut String,
    distances: &[f64],
    fluxes: &[f64],
    d_min: f64,
    d_range: f64,
    f_min_log: f64,
    f_range_log: f64,
    margin: f64,
    plot_w: f64,
    plot_h: f64,
    color: &str,
) {
    let mut points = String::new();
    for (i, (d, f)) in distances.iter().zip(fluxes.iter()).enumerate() {
        if *f <= 0.0 {
            continue;
        }
        let x = margin + (d - d_min) / d_range * plot_w;
        let y = margin + plot_h - ((f.log10() - f_min_log) / f_range_log * plot_h);
        if i == 0 {
            write!(points, "{x:.1},{y:.1}").unwrap();
        } else {
            write!(points, " {x:.1},{y:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>"#,
    )
    .unwrap();
}

/// Generate an SVG diagram showing the angular separation vs. distance
/// relationship, highlighting the 42'–69.6' search window.
pub fn separation_vs_distance_plot(
    distances: &[f64],
    separations: &[f64],
    width: u32,
    height: u32,
) -> String {
    assert_eq!(distances.len(), separations.len());
    if distances.is_empty() {
        return String::new();
    }

    let margin = 70.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let d_min = distances.first().copied().unwrap_or(200.0);
    let d_max = distances.last().copied().unwrap_or(1000.0);
    let d_range = d_max - d_min;

    let s_min = 0.0_f64;
    let s_max = separations.iter().copied().fold(0.0_f64, f64::max) * 1.2;
    let s_range = s_max - s_min;

    let mut svg = String::with_capacity(3000);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Angular Separation over 23-yr Baseline</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();
    writeln!(
        svg,
        r##"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="#f8f8f8" stroke="black" stroke-width="1"/>"##,
    )
    .unwrap();

    // Search window horizontal lines at 42' and 69.6'
    for &(sep_val, label) in &[(42.0, "42' (700 AU)"), (69.6, "69.6' (500 AU)")] {
        let y = margin + plot_h - (sep_val - s_min) / s_range * plot_h;
        writeln!(
            svg,
            r#"<line x1="{margin}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="goldenrod" stroke-width="1" stroke-dasharray="4,4"/>"#,
            margin + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{}" font-size="8" font-family="sans-serif" fill="goldenrod">{label}</text>"#,
            margin + plot_w + 3.0,
            y + 3.0,
        )
        .unwrap();
    }

    // Shaded region between the two lines
    let y_42 = margin + plot_h - (42.0 - s_min) / s_range * plot_h;
    let y_69 = margin + plot_h - (69.6 - s_min) / s_range * plot_h;
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{y_69:.1}" width="{plot_w}" height="{:.1}" fill="gold" fill-opacity="0.15"/>"#,
        y_42 - y_69,
    )
    .unwrap();

    // Curve
    let mut points = String::new();
    for (i, (d, s)) in distances.iter().zip(separations.iter()).enumerate() {
        let x = margin + (d - d_min) / d_range * plot_w;
        let y = margin + plot_h - (s - s_min) / s_range * plot_h;
        if i == 0 {
            write!(points, "{x:.1},{y:.1}").unwrap();
        } else {
            write!(points, " {x:.1},{y:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<polyline points="{points}" fill="none" stroke="steelblue" stroke-width="2"/>"#,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Heliocentric Distance (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 40.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Angular Separation (arcmin)</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // X-axis ticks
    for &tick in &[200.0, 400.0, 600.0, 800.0, 1000.0] {
        if tick >= d_min && tick <= d_max {
            let x = margin + (tick - d_min) / d_range * plot_w;
            writeln!(
                svg,
                r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="8" font-family="sans-serif">{tick:.0}</text>"#,
                margin + plot_h + 15.0,
            )
            .unwrap();
        }
    }

    // Y-axis ticks
    let y_step = if s_max > 100.0 { 20.0 } else { 10.0 };
    let mut tick = 0.0;
    while tick <= s_max {
        let y = margin + plot_h - (tick - s_min) / s_range * plot_h;
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="8" font-family="sans-serif">{tick:.0}'</text>"#,
            margin - 5.0,
        )
        .unwrap();
        tick += y_step;
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::coords::candidate_pair::heliocentric_separation_arcmin;

    #[test]
    fn flux_plot_produces_valid_svg() {
        let distances: Vec<f64> = (200..=1000).step_by(100).map(|d| d as f64).collect();
        let flux_60: Vec<f64> = distances.iter().map(|d| 1.0 / (d * d) * 1e5).collect();
        let flux_90: Vec<f64> = distances.iter().map(|d| 1.5 / (d * d) * 1e5).collect();
        let svg = flux_vs_distance_plot(&distances, &flux_60, &flux_90, 700, 450);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("IRAS"));
        assert!(svg.contains("AKARI"));
        assert!(svg.contains("Search Window"));
    }

    #[test]
    fn separation_plot_produces_valid_svg() {
        let distances: Vec<f64> = (200..=1000).step_by(50).map(|d| d as f64).collect();
        let seps: Vec<f64> = distances
            .iter()
            .map(|d| heliocentric_separation_arcmin(*d, *d, 0.0, 23.0))
            .collect();
        let svg = separation_vs_distance_plot(&distances, &seps, 700, 450);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("42'"));
        assert!(svg.contains("69.6'"));
    }
}
