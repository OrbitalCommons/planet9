//! SVG plot generation for DES Planet Nine analysis figures.
//!
//! Produces visualizations of color model recovery rates and
//! detection efficiency as a function of magnitude.

use std::fmt::Write;

use crate::color_models;
use crate::survey_model::DesSurvey;

const GRIDLINE_COLOR: &str = "silver";
const CURVE_COLOR: &str = "steelblue";

/// Generate a bar chart comparing recovery rates across color models.
///
/// Each bar represents a different surface color/albedo assumption
/// from Table 1, with the recovery rate on the y-axis.
pub fn color_model_comparison_plot(width: u32, height: u32) -> String {
    let models = color_models::all_models();
    let margin_left = 70.0_f64;
    let margin_right = 20.0_f64;
    let margin_top = 40.0_f64;
    let margin_bottom = 80.0_f64;
    let plot_w = width as f64 - margin_left - margin_right;
    let plot_h = height as f64 - margin_top - margin_bottom;

    let n = models.len();
    let bar_width = plot_w / (n as f64 * 1.5 + 0.5);
    let bar_spacing = bar_width * 0.5;

    let colors = [
        "steelblue",
        "darkorange",
        "indianred",
        "cadetblue",
        "olivedrab",
    ];

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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">DES Recovery Rate by Color Model</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Y-axis gridlines and labels
    for &frac in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = margin_top + plot_h * (1.0 - frac);
        writeln!(
            svg,
            r#"<line x1="{margin_left}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="{GRIDLINE_COLOR}" stroke-width="0.5"/>"#,
            margin_left + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{}" text-anchor="end" font-size="10" font-family="sans-serif" dominant-baseline="middle">{:.0}%</text>"#,
            margin_left - 5.0,
            y,
            frac * 100.0,
        )
        .unwrap();
    }

    // Y-axis label
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Recovery Rate</text>"#,
        margin_top + plot_h / 2.0,
        margin_top + plot_h / 2.0,
    )
    .unwrap();

    // Bars
    for (i, model) in models.iter().enumerate() {
        let x = margin_left + bar_spacing + i as f64 * (bar_width + bar_spacing);
        let bar_h = plot_h * model.recovery_rate;
        let y = margin_top + plot_h - bar_h;
        let color = colors[i % colors.len()];

        writeln!(
            svg,
            r#"<rect x="{x:.1}" y="{y:.1}" width="{bar_width:.1}" height="{bar_h:.1}" fill="{color}" stroke="black" stroke-width="0.5"/>"#,
        )
        .unwrap();

        // Value label
        writeln!(
            svg,
            r#"<text x="{}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{:.0}%</text>"#,
            x + bar_width / 2.0,
            y - 5.0,
            model.recovery_rate * 100.0,
        )
        .unwrap();

        // X-axis label (rotated)
        writeln!(
            svg,
            r#"<text x="{}" y="{}" text-anchor="end" font-size="9" font-family="sans-serif" transform="rotate(-45, {}, {})">{}</text>"#,
            x + bar_width / 2.0,
            margin_top + plot_h + 15.0,
            x + bar_width / 2.0,
            margin_top + plot_h + 15.0,
            model.name,
        )
        .unwrap();
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a detection efficiency curve as a function of apparent magnitude.
///
/// Shows the DES completeness (logistic function) from bright to faint,
/// with the survey depth limit marked.
pub fn recovery_vs_magnitude_plot(width: u32, height: u32) -> String {
    let survey = DesSurvey::default();
    let margin_left = 60.0_f64;
    let margin_right = 20.0_f64;
    let margin_top = 40.0_f64;
    let margin_bottom = 50.0_f64;
    let plot_w = width as f64 - margin_left - margin_right;
    let plot_h = height as f64 - margin_top - margin_bottom;

    let mag_min = 18.0_f64;
    let mag_max = 27.0_f64;
    let mag_range = mag_max - mag_min;

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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">DES Detection Efficiency vs. Magnitude</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Y-axis gridlines and labels
    for &frac in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = margin_top + plot_h * (1.0 - frac);
        writeln!(
            svg,
            r#"<line x1="{margin_left}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="{GRIDLINE_COLOR}" stroke-width="0.5"/>"#,
            margin_left + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{}" text-anchor="end" font-size="10" font-family="sans-serif" dominant-baseline="middle">{:.0}%</text>"#,
            margin_left - 5.0,
            y,
            frac * 100.0,
        )
        .unwrap();
    }

    // X-axis labels
    for mag in (18..=27).step_by(1) {
        let x = margin_left + plot_w * (mag as f64 - mag_min) / mag_range;
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="10" font-family="sans-serif">{mag}</text>"#,
            margin_top + plot_h + 18.0,
        )
        .unwrap();
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Apparent g-band Magnitude</text>"#,
        margin_left + plot_w / 2.0,
        margin_top + plot_h + 38.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Completeness</text>"#,
        margin_top + plot_h / 2.0,
        margin_top + plot_h / 2.0,
    )
    .unwrap();

    // Efficiency curve
    let n_points = 200;
    let mut path = String::with_capacity(2000);
    for i in 0..=n_points {
        let mag = mag_min + mag_range * (i as f64 / n_points as f64);
        let eff = survey.completeness(mag);
        let x = margin_left + plot_w * (mag - mag_min) / mag_range;
        let y = margin_top + plot_h * (1.0 - eff);
        if i == 0 {
            write!(path, "M{x:.1},{y:.1}").unwrap();
        } else {
            write!(path, " L{x:.1},{y:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<path d="{path}" fill="none" stroke="{CURVE_COLOR}" stroke-width="2"/>"#,
    )
    .unwrap();

    // Depth limit vertical line
    let depth_x = margin_left + plot_w * (survey.depth_g - mag_min) / mag_range;
    writeln!(
        svg,
        r#"<line x1="{depth_x:.1}" y1="{margin_top}" x2="{depth_x:.1}" y2="{}" stroke="red" stroke-width="1" stroke-dasharray="4,4"/>"#,
        margin_top + plot_h,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="start" font-size="9" font-family="sans-serif" fill="red">g = {:.1}</text>"#,
        depth_x + 4.0,
        margin_top + 15.0,
        survey.depth_g,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn color_model_plot_is_valid_svg() {
        let svg = color_model_comparison_plot(600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("DES Recovery Rate"));
    }

    #[test]
    fn color_model_plot_contains_all_models() {
        let svg = color_model_comparison_plot(600, 400);
        for model in color_models::all_models() {
            assert!(
                svg.contains(&model.name),
                "Plot missing model: {}",
                model.name
            );
        }
    }

    #[test]
    fn magnitude_plot_is_valid_svg() {
        let svg = recovery_vs_magnitude_plot(600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Detection Efficiency"));
    }

    #[test]
    fn magnitude_plot_contains_depth_marker() {
        let svg = recovery_vs_magnitude_plot(600, 400);
        assert!(svg.contains("24.1"));
    }
}
