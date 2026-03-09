//! SVG plot generation for resonance search figures.
//!
//! Reproduces key figures from Bailey+ (2018):
//! - Resonance histogram (Figure 4)
//! - a₉ distribution comparison (Figure 6)

use std::fmt::Write;

use crate::probability_analysis::{A9Distribution, DistributionComparison};
use crate::resonance_catalog::Resonance;

/// Generate a histogram of resonance occupation.
///
/// Shows how many particles occupy each resonance, colored by type
/// (N/1 = red, N/2 = orange, high-order = blue).
pub fn resonance_histogram_plot(
    census: &[(Resonance, usize)],
    e_p9: f64,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let max_count = census.iter().map(|(_, c)| *c).max().unwrap_or(1) as f64;
    let n_bars = census.len();

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Resonance Occupation (e₉ = {e_p9:.1})</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    if n_bars > 0 {
        let bar_w = (plot_w / n_bars as f64).min(30.0);

        for (i, (res, count)) in census.iter().enumerate() {
            let x = margin + (i as f64 + 0.5) * (plot_w / n_bars as f64) - bar_w / 2.0;
            let h = (*count as f64 / max_count) * plot_h;
            let y = margin + plot_h - h;

            let color = if res.is_n_over_1() {
                "indianred"
            } else if res.is_n_over_2() {
                "orange"
            } else {
                "steelblue"
            };

            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{bar_w:.1}" height="{h:.1}" fill="{color}" stroke="black" stroke-width="0.3"/>"#,
            )
            .unwrap();

            // Label
            let label_x = margin + (i as f64 + 0.5) * (plot_w / n_bars as f64);
            writeln!(
                svg,
                r#"<text x="{label_x:.1}" y="{}" text-anchor="middle" font-size="7" font-family="sans-serif" transform="rotate(-45, {label_x:.1}, {})">{res}</text>"#,
                margin + plot_h + 12.0,
                margin + plot_h + 12.0,
            )
            .unwrap();
        }
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Resonance p:q</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 40.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Count</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Legend
    let lx = margin + plot_w - 130.0;
    let ly = margin + 10.0;
    for (i, (color, label)) in [
        ("indianred", "N:1"),
        ("orange", "N:2"),
        ("steelblue", "High-order"),
    ]
    .iter()
    .enumerate()
    {
        let y = ly + i as f64 * 16.0;
        writeln!(
            svg,
            r#"<rect x="{lx}" y="{y}" width="12" height="12" fill="{color}" stroke="black" stroke-width="0.3"/>"#,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif">{label}</text>"#,
            lx + 16.0,
            y + 10.0,
        )
        .unwrap();
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a comparison plot of a₉ distributions (Figure 6).
///
/// Shows Farey F5 (peaked) vs extended catalog (broad plateau).
pub fn a9_comparison_plot(comparison: &DistributionComparison, width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let max_density = comparison
        .f5_distribution
        .density
        .iter()
        .chain(comparison.ext_distribution.density.iter())
        .copied()
        .fold(0.0_f64, f64::max);

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">a₉ Distribution: F5 vs Extended Resonances</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Plot both distributions as lines
    if max_density > 0.0 {
        draw_distribution_line(
            &mut svg,
            &comparison.f5_distribution,
            margin,
            plot_w,
            plot_h,
            max_density,
            "indianred",
        );
        draw_distribution_line(
            &mut svg,
            &comparison.ext_distribution,
            margin,
            plot_w,
            plot_h,
            max_density,
            "steelblue",
        );
    }

    // X-axis labels
    let a_range = (200.0, 1200.0);
    for a_val in [200, 400, 600, 800, 1000, 1200] {
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
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Planet Nine Semi-major Axis a₉ (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Probability Density</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Legend
    let lx = margin + plot_w - 180.0;
    let ly = margin + 10.0;
    writeln!(
        svg,
        r#"<line x1="{lx}" y1="{}" x2="{}" y2="{}" stroke="indianred" stroke-width="2"/>"#,
        ly + 6.0,
        lx + 20.0,
        ly + 6.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif">Farey F5 (peak/mean={:.1})</text>"#,
        lx + 25.0,
        ly + 10.0,
        comparison.f5_peak_to_mean,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<line x1="{lx}" y1="{}" x2="{}" y2="{}" stroke="steelblue" stroke-width="2"/>"#,
        ly + 22.0,
        lx + 20.0,
        ly + 22.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif">Extended (peak/mean={:.1})</text>"#,
        lx + 25.0,
        ly + 26.0,
        comparison.ext_peak_to_mean,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

fn draw_distribution_line(
    svg: &mut String,
    dist: &A9Distribution,
    margin: f64,
    plot_w: f64,
    plot_h: f64,
    max_density: f64,
    color: &str,
) {
    let a_range = (200.0, 1200.0);
    let mut points = String::new();

    for (i, (&bin, &d)) in dist.bins.iter().zip(dist.density.iter()).enumerate() {
        let x = margin + (bin - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        let y = margin + plot_h - (d / max_density) * plot_h;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::probability_analysis::{compare_distributions, OBSERVED_KBO_AXES};

    #[test]
    fn test_resonance_histogram_plot() {
        let census = vec![
            (Resonance::new(2, 1), 15),
            (Resonance::new(3, 2), 8),
            (Resonance::new(5, 3), 3),
        ];

        let svg = resonance_histogram_plot(&census, 0.3, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Resonance Occupation"));
    }

    #[test]
    fn test_a9_comparison_plot() {
        let comparison = compare_distributions(OBSERVED_KBO_AXES);
        let svg = a9_comparison_plot(&comparison, 700, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Farey F5"));
        assert!(svg.contains("Extended"));
    }
}
