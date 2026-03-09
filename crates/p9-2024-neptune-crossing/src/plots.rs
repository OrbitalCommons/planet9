//! SVG plot generation for Neptune-crossing TNO analysis figures.
//!
//! Reproduces key figures from the paper:
//! - CDF comparison of xi values for P9-inclusive vs P9-free vs uniform
//! - Perihelion distance distribution comparison

use std::fmt::Write;

use crate::simulation::SimulationResult;

/// Generate a CDF comparison plot of xi values.
///
/// Plots the empirical CDFs of the CDF-transformed perihelion values (xi)
/// for the P9-inclusive model, the P9-free model, and the expected uniform
/// distribution.
pub fn cdf_comparison_plot(xi_p9: &[f64], xi_null: &[f64], width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">CDF of Perihelion Test Statistic</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Uniform diagonal (expected under correct model)
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{}" x2="{}" y2="{margin}" stroke="gray" stroke-width="1" stroke-dasharray="6,3"/>"#,
        margin + plot_h,
        margin + plot_w,
    )
    .unwrap();

    // Draw empirical CDF for P9 model (blue)
    draw_ecdf(&mut svg, xi_p9, margin, plot_w, plot_h, "steelblue", 2.0);

    // Draw empirical CDF for null model (red)
    draw_ecdf(&mut svg, xi_null, margin, plot_w, plot_h, "crimson", 2.0);

    // Legend
    let lx = margin + plot_w - 130.0;
    let ly = margin + 25.0;
    writeln!(
        svg,
        r#"<line x1="{lx}" y1="{ly}" x2="{}" y2="{ly}" stroke="steelblue" stroke-width="2"/>"#,
        lx + 20.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">P9-inclusive</text>"#,
        lx + 25.0,
        ly + 4.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<line x1="{lx}" y1="{}" x2="{}" y2="{}" stroke="crimson" stroke-width="2"/>"#,
        ly + 18.0,
        lx + 20.0,
        ly + 18.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">P9-free</text>"#,
        lx + 25.0,
        ly + 22.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<line x1="{lx}" y1="{}" x2="{}" y2="{}" stroke="gray" stroke-width="1" stroke-dasharray="6,3"/>"#,
        ly + 36.0,
        lx + 20.0,
        ly + 36.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">Uniform</text>"#,
        lx + 25.0,
        ly + 40.0,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">xi</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">Cumulative Fraction</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a perihelion distance distribution comparison plot.
///
/// Shows the distribution of perihelion distances for P9-inclusive vs P9-free
/// simulation results alongside observed TNO perihelia.
pub fn perihelion_distribution_plot(
    p9_result: &SimulationResult,
    null_result: &SimulationResult,
    observed_q: &[f64],
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let q_min = 10.0_f64;
    let q_max = 32.0_f64;
    let n_bins = 11;
    let bin_width = (q_max - q_min) / n_bins as f64;

    let mut svg = String::with_capacity(5000);
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Perihelion Distribution (q &lt; 30 AU)</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Bin the data
    let p9_bins = bin_data(&p9_result.selected_perihelia, q_min, q_max, n_bins);
    let null_bins = bin_data(&null_result.selected_perihelia, q_min, q_max, n_bins);
    let obs_bins = bin_data(observed_q, q_min, q_max, n_bins);

    // Normalize to peak = 1 for display
    let p9_max = *p9_bins.iter().max().unwrap_or(&1) as f64;
    let null_max = *null_bins.iter().max().unwrap_or(&1) as f64;
    let obs_max = *obs_bins.iter().max().unwrap_or(&1) as f64;
    let global_max = p9_max.max(null_max).max(obs_max).max(1.0);

    let bar_group_width = plot_w / n_bins as f64;
    let bar_width = bar_group_width * 0.25;

    for i in 0..n_bins {
        let x_base = margin + i as f64 * bar_group_width;

        // P9 bars (blue)
        let h_p9 = (p9_bins[i] as f64 / global_max) * plot_h;
        writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="steelblue" opacity="0.7"/>"#,
            x_base + 2.0,
            margin + plot_h - h_p9,
            bar_width,
            h_p9,
        )
        .unwrap();

        // Null bars (red)
        let h_null = (null_bins[i] as f64 / global_max) * plot_h;
        writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="crimson" opacity="0.7"/>"#,
            x_base + 2.0 + bar_width,
            margin + plot_h - h_null,
            bar_width,
            h_null,
        )
        .unwrap();

        // Observed bars (green)
        let h_obs = (obs_bins[i] as f64 / global_max) * plot_h;
        writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="forestgreen" opacity="0.7"/>"#,
            x_base + 2.0 + 2.0 * bar_width,
            margin + plot_h - h_obs,
            bar_width,
            h_obs,
        )
        .unwrap();

        // Bin label
        let bin_center = q_min + (i as f64 + 0.5) * bin_width;
        writeln!(
            svg,
            r#"<text x="{:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{:.0}</text>"#,
            x_base + bar_group_width / 2.0,
            margin + plot_h + 15.0,
            bin_center,
        )
        .unwrap();
    }

    // Legend
    let lx = margin + 10.0;
    let ly = margin + 15.0;
    writeln!(
        svg,
        r#"<rect x="{lx}" y="{}" width="12" height="12" fill="steelblue" opacity="0.7"/>"#,
        ly - 10.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{ly}" font-size="10" font-family="sans-serif">P9-inclusive</text>"#,
        lx + 16.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect x="{lx}" y="{}" width="12" height="12" fill="crimson" opacity="0.7"/>"#,
        ly + 6.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">P9-free</text>"#,
        lx + 16.0,
        ly + 16.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect x="{lx}" y="{}" width="12" height="12" fill="forestgreen" opacity="0.7"/>"#,
        ly + 22.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">Observed</text>"#,
        lx + 16.0,
        ly + 32.0,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Perihelion q (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">Count</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

fn draw_ecdf(
    svg: &mut String,
    values: &[f64],
    margin: f64,
    plot_w: f64,
    plot_h: f64,
    color: &str,
    stroke_width: f64,
) {
    if values.is_empty() {
        return;
    }

    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len() as f64;

    let mut points = Vec::with_capacity(sorted.len() + 2);
    points.push(format!("{:.1},{:.1}", margin, margin + plot_h,));

    for (i, &val) in sorted.iter().enumerate() {
        let x = margin + val.clamp(0.0, 1.0) * plot_w;
        let y = margin + plot_h - ((i + 1) as f64 / n) * plot_h;
        points.push(format!("{x:.1},{y:.1}"));
    }

    writeln!(
        svg,
        r#"<polyline points="{}" fill="none" stroke="{color}" stroke-width="{stroke_width}"/>"#,
        points.join(" "),
    )
    .unwrap();
}

fn bin_data(values: &[f64], min: f64, max: f64, n_bins: usize) -> Vec<usize> {
    let mut bins = vec![0usize; n_bins];
    let bin_width = (max - min) / n_bins as f64;

    for &v in values {
        if v >= min && v < max {
            let idx = ((v - min) / bin_width) as usize;
            let idx = idx.min(n_bins - 1);
            bins[idx] += 1;
        }
    }
    bins
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::quick_test_simulation;

    #[test]
    fn test_cdf_plot_valid_svg() {
        let xi_p9 = vec![0.1, 0.3, 0.45, 0.5, 0.6, 0.8, 0.95];
        let xi_null = vec![0.02, 0.05, 0.12, 0.2, 0.35, 0.5, 0.7];
        let svg = cdf_comparison_plot(&xi_p9, &xi_null, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("steelblue"));
        assert!(svg.contains("crimson"));
    }

    #[test]
    fn test_perihelion_distribution_plot_valid_svg() {
        let p9_result = quick_test_simulation(true);
        let null_result = quick_test_simulation(false);
        let observed_q: Vec<f64> = crate::observed_tnos::observed_sample()
            .iter()
            .map(|t| t.q)
            .collect();

        let svg = perihelion_distribution_plot(&p9_result, &null_result, &observed_q, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Perihelion Distribution"));
    }

    #[test]
    fn test_cdf_plot_empty_data() {
        let svg = cdf_comparison_plot(&[], &[], 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_bin_data() {
        let values = vec![11.0, 15.0, 15.5, 20.0, 25.0, 29.0];
        let bins = bin_data(&values, 10.0, 30.0, 4);
        assert_eq!(bins.len(), 4);
        let total: usize = bins.iter().sum();
        assert_eq!(total, 6);
    }
}
