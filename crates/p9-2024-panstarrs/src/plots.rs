//! SVG plot generation for paper figures.
//!
//! Generates SVG visualizations reproducing the key figures from
//! Brown, Holman & Batygin (2024):
//! - Combined exclusion stacked bar chart
//! - Remaining parameter space schematic

use std::fmt::Write;

use crate::combined_exclusion::CombinedExclusion;

/// Generate a stacked bar chart showing exclusion fraction by survey.
///
/// Produces an SVG with three stacked segments (ZTF, DES, PS1) showing
/// each survey's contribution to the total 78% exclusion.
pub fn combined_exclusion_plot(exclusion: &CombinedExclusion, width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let bar_width = 120.0_f64;
    let plot_h = height as f64 - 2.0 * margin;
    let bar_x = (width as f64 - bar_width) / 2.0;

    let mut svg = String::with_capacity(3000);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#
    )
    .unwrap();

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Combined Exclusion by Survey</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Y-axis: 0% at bottom, 100% at top
    let y_for_frac = |frac: f64| -> f64 { margin + plot_h * (1.0 - frac) };

    // Stacked segments (bottom to top): ZTF, DES unique, PS1 unique
    let segments: Vec<(&str, f64, f64, &str)> = vec![
        ("ZTF (56.4%)", 0.0, exclusion.ztf_frac, "steelblue"),
        (
            "DES (5.0%)",
            exclusion.ztf_frac,
            exclusion.ztf_frac + exclusion.des_unique,
            "indianred",
        ),
        (
            "PS1 (17.1%)",
            exclusion.ztf_frac + exclusion.des_unique,
            exclusion.combined,
            "seagreen",
        ),
    ];

    for (label, frac_start, frac_end, color) in &segments {
        let y_top = y_for_frac(*frac_end);
        let y_bot = y_for_frac(*frac_start);
        let seg_h = y_bot - y_top;

        writeln!(
            svg,
            r#"<rect x="{bar_x}" y="{y_top:.1}" width="{bar_width}" height="{seg_h:.1}" fill="{color}" stroke="white" stroke-width="1"/>"#,
        )
        .unwrap();

        // Label to the right of the bar
        let label_y = y_top + seg_h / 2.0 + 4.0;
        writeln!(
            svg,
            r#"<text x="{}" y="{label_y:.1}" font-size="11" font-family="sans-serif">{label}</text>"#,
            bar_x + bar_width + 8.0,
        )
        .unwrap();
    }

    // Remaining space
    let remaining_y_top = y_for_frac(1.0);
    let remaining_y_bot = y_for_frac(exclusion.combined);
    let remaining_h = remaining_y_bot - remaining_y_top;
    let gray = "silver";
    writeln!(
        svg,
        r#"<rect x="{bar_x}" y="{remaining_y_top:.1}" width="{bar_width}" height="{remaining_h:.1}" fill="{gray}" stroke="white" stroke-width="1"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{:.1}" font-size="11" font-family="sans-serif">Remaining (22%)</text>"#,
        bar_x + bar_width + 8.0,
        remaining_y_top + remaining_h / 2.0 + 4.0,
    )
    .unwrap();

    // Y-axis ticks
    for pct in [0, 25, 50, 75, 100] {
        let frac = pct as f64 / 100.0;
        let y = y_for_frac(frac);
        writeln!(
            svg,
            r#"<line x1="{}" x2="{bar_x}" y1="{y:.1}" y2="{y:.1}" stroke="black" stroke-width="0.5"/>"#,
            bar_x - 5.0,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{:.1}" text-anchor="end" font-size="10" font-family="sans-serif">{pct}%</text>"#,
            bar_x - 8.0,
            y + 3.5,
        )
        .unwrap();
    }

    // Axis line
    writeln!(
        svg,
        r#"<line x1="{bar_x}" x2="{bar_x}" y1="{margin}" y2="{:.1}" stroke="black" stroke-width="1"/>"#,
        margin + plot_h,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a schematic of the remaining viable parameter space.
///
/// Shows a semi-major axis vs V-magnitude scatter region with the excluded
/// zone shaded and the remaining 22% highlighted.
pub fn remaining_sky_plot(width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(3000);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#
    )
    .unwrap();

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Remaining Parameter Space</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot area
    let bg = "gainsboro";
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="{bg}" stroke="black"/>"#,
    )
    .unwrap();

    // Excluded region (brighter and closer objects already surveyed)
    // a: 200-800 AU mapped to x; V: 19-25 mapped to y (bright at top)
    let a_min = 200.0_f64;
    let a_max = 800.0_f64;
    let v_min = 19.0_f64;
    let v_max = 25.0_f64;

    let x_for_a = |a: f64| -> f64 { margin + (a - a_min) / (a_max - a_min) * plot_w };
    let y_for_v = |v: f64| -> f64 { margin + (v - v_min) / (v_max - v_min) * plot_h };

    // ZTF excluded region: bright objects (V < 20.5) across most of sky
    let ztf_color = "steelblue";
    writeln!(
        svg,
        r#"<rect x="{}" y="{}" width="{}" height="{:.1}" fill="{ztf_color}" opacity="0.4"/>"#,
        x_for_a(a_min),
        y_for_v(v_min),
        x_for_a(a_max) - x_for_a(a_min),
        y_for_v(20.5) - y_for_v(v_min),
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="10" font-family="sans-serif" fill="{ztf_color}">ZTF excluded</text>"#,
        x_for_a(250.0),
        y_for_v(19.8),
    )
    .unwrap();

    // PS1 extends to fainter magnitudes
    let ps1_color = "seagreen";
    writeln!(
        svg,
        r#"<rect x="{}" y="{:.1}" width="{}" height="{:.1}" fill="{ps1_color}" opacity="0.3"/>"#,
        x_for_a(a_min),
        y_for_v(20.5),
        x_for_a(600.0) - x_for_a(a_min),
        y_for_v(21.5) - y_for_v(20.5),
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="10" font-family="sans-serif" fill="{ps1_color}">PS1 excluded</text>"#,
        x_for_a(250.0),
        y_for_v(21.1),
    )
    .unwrap();

    // Remaining region marker (faint + distant)
    let highlight = "crimson";
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="12" font-family="sans-serif" fill="{highlight}" font-weight="bold">Remaining 22%</text>"#,
        x_for_a(450.0),
        y_for_v(23.0),
    )
    .unwrap();

    // Updated parameter estimate marker
    let est_x = x_for_a(500.0);
    let est_y = y_for_v(22.0);
    writeln!(
        svg,
        r#"<circle cx="{est_x:.1}" cy="{est_y:.1}" r="6" fill="none" stroke="{highlight}" stroke-width="2"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="9" font-family="sans-serif" fill="{highlight}">a=500 AU, V=22.0</text>"#,
        est_x + 10.0,
        est_y + 4.0,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Semi-major axis (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 10.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">V magnitude</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    // X-axis ticks
    for a in [200, 400, 600, 800] {
        let x = x_for_a(a as f64);
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{:.1}" text-anchor="middle" font-size="10" font-family="sans-serif">{a}</text>"#,
            margin + plot_h + 15.0,
        )
        .unwrap();
    }

    // Y-axis ticks
    for v in [19, 20, 21, 22, 23, 24, 25] {
        let y = y_for_v(v as f64);
        writeln!(
            svg,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="end" font-size="10" font-family="sans-serif">{v}</text>"#,
            margin - 5.0,
            y + 3.5,
        )
        .unwrap();
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exclusion_plot_valid_svg() {
        let ex = CombinedExclusion::paper_values();
        let svg = combined_exclusion_plot(&ex, 500, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("rect"));
    }

    #[test]
    fn exclusion_plot_contains_survey_labels() {
        let ex = CombinedExclusion::paper_values();
        let svg = combined_exclusion_plot(&ex, 500, 400);
        assert!(svg.contains("ZTF"));
        assert!(svg.contains("DES"));
        assert!(svg.contains("PS1"));
        assert!(svg.contains("Remaining"));
    }

    #[test]
    fn remaining_sky_plot_valid_svg() {
        let svg = remaining_sky_plot(600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("circle"));
    }

    #[test]
    fn remaining_sky_plot_contains_labels() {
        let svg = remaining_sky_plot(600, 400);
        assert!(svg.contains("Remaining 22%"));
        assert!(svg.contains("a=500 AU"));
        assert!(svg.contains("Semi-major axis"));
        assert!(svg.contains("V magnitude"));
    }
}
