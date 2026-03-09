//! SVG plot generation for the review paper figures.
//!
//! Focuses on:
//! - Parameter space comparison (2016 vs 2019)
//! - Brightness curves with survey depth overlays

use std::fmt::Write;

use crate::detection_prospects::{
    apparent_magnitude, brightness_table, survey_limits, BrightnessEstimate, P9PhysicalProperties,
};
use crate::revised_parameters::{original_2016, revised_2019, P9ParameterSet};

/// Generate a parameter comparison chart (2016 vs 2019).
pub fn parameter_comparison_plot(width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let orig = original_2016();
    let rev = revised_2019();

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Planet Nine Parameters: 2016 vs 2019</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot in (a, e) space
    let a_range = (200.0, 1000.0);
    let e_range = (0.0, 0.8);

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // 2016 point
    let x_16 = margin + (orig.a.best - a_range.0) / (a_range.1 - a_range.0) * plot_w;
    let y_16 = margin + plot_h - (orig.e.best - e_range.0) / (e_range.1 - e_range.0) * plot_h;
    writeln!(
        svg,
        r#"<circle cx="{x_16:.1}" cy="{y_16:.1}" r="8" fill="none" stroke="indianred" stroke-width="2"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif" fill="indianred">2016</text>"#,
        x_16 + 12.0,
        y_16 + 4.0,
    )
    .unwrap();

    // 2019 range (rectangle)
    let x_min = margin + (rev.a.min - a_range.0) / (a_range.1 - a_range.0) * plot_w;
    let x_max = margin + (rev.a.max - a_range.0) / (a_range.1 - a_range.0) * plot_w;
    let y_min = margin + plot_h - (rev.e.max - e_range.0) / (e_range.1 - e_range.0) * plot_h;
    let y_max = margin + plot_h - (rev.e.min - e_range.0) / (e_range.1 - e_range.0) * plot_h;
    writeln!(
        svg,
        r#"<rect x="{x_min:.1}" y="{y_min:.1}" width="{:.1}" height="{:.1}" fill="steelblue" fill-opacity="0.2" stroke="steelblue" stroke-width="1.5"/>"#,
        x_max - x_min,
        y_max - y_min,
    )
    .unwrap();

    // 2019 best-fit point
    let x_19 = margin + (rev.a.best - a_range.0) / (a_range.1 - a_range.0) * plot_w;
    let y_19 = margin + plot_h - (rev.e.best - e_range.0) / (e_range.1 - e_range.0) * plot_h;
    writeln!(
        svg,
        r#"<circle cx="{x_19:.1}" cy="{y_19:.1}" r="6" fill="steelblue" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif" fill="steelblue">2019</text>"#,
        x_19 + 10.0,
        y_19 + 4.0,
    )
    .unwrap();

    // Axis labels
    for a_val in [200, 400, 600, 800, 1000] {
        let x = margin + (a_val as f64 - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{a_val}</text>"#,
            margin + plot_h + 15.0,
        )
        .unwrap();
    }
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Semi-major Axis a₉ (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    for e_val in [0.0, 0.2, 0.4, 0.6, 0.8] {
        let y = margin + plot_h - (e_val as f64 - e_range.0) / (e_range.1 - e_range.0) * plot_h;
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="9" font-family="sans-serif">{e_val:.1}</text>"#,
            margin - 5.0,
        )
        .unwrap();
    }
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Eccentricity e₉</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a brightness curve plot.
pub fn brightness_plot(width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let r_range = (100.0, 1200.0);
    let v_range = (16.0, 28.0);

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Planet Nine Apparent Brightness</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Draw brightness curves for 5 ME and 10 ME (optimistic)
    let configs = [
        (
            P9PhysicalProperties::five_me_optimistic(),
            "steelblue",
            "5 ME (bright)",
        ),
        (
            P9PhysicalProperties::ten_me_conservative(),
            "indianred",
            "10 ME (faint)",
        ),
    ];

    for (physical, color, _label) in &configs {
        let mut points = String::new();
        for i in 0..100 {
            let r = r_range.0 + (r_range.1 - r_range.0) * i as f64 / 99.0;
            let v = apparent_magnitude(r, physical);
            let x = margin + (r - r_range.0) / (r_range.1 - r_range.0) * plot_w;
            let y = margin + (v - v_range.0) / (v_range.1 - v_range.0) * plot_h;
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

    // Survey depth lines
    for limit in survey_limits() {
        let y = margin + (limit.v_limit - v_range.0) / (v_range.1 - v_range.0) * plot_h;
        if y > margin && y < margin + plot_h {
            writeln!(
                svg,
                r#"<line x1="{margin}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="gray" stroke-width="1" stroke-dasharray="4,4"/>"#,
                margin + plot_w,
            )
            .unwrap();
            writeln!(
                svg,
                r#"<text x="{}" y="{}" font-size="8" font-family="sans-serif" fill="gray">{}</text>"#,
                margin + plot_w - 2.0,
                y - 3.0,
                limit.name,
            )
            .unwrap();
        }
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Heliocentric Distance (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">V Magnitude</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parameter_comparison_plot() {
        let svg = parameter_comparison_plot(600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("2016"));
        assert!(svg.contains("2019"));
    }

    #[test]
    fn test_brightness_plot() {
        let svg = brightness_plot(700, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Brightness"));
    }
}
