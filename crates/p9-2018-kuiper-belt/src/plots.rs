//! SVG plot generation for Kuiper Belt generation figures.
//!
//! Reproduces Figure 1/2 from Khain+ 2018: side-by-side perihelion
//! distributions for narrow vs broad initial conditions, colored by
//! alignment with Planet Nine (classified from the Δϖ(t) snapshot series,
//! not a single snapshot).

use std::fmt::Write;

use crate::population_analysis::{
    classify_alignment_series, perihelion_histograms_by_alignment, population_statistics_series,
    Alignment, PopulationResult,
};
use crate::simulation::SimulationResult;

/// Generate a perihelion distance histogram SVG colored by alignment.
///
/// Reproduces Figure 1 of Khain+ 2018: distribution of final q values
/// with aligned (red) and anti-aligned (blue) populations.
pub fn perihelion_distribution_plot(
    result: &SimulationResult,
    title: &str,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let q_range = (0.0, 350.0);
    let n_bins = 35;
    let (centers, aligned, anti_aligned) =
        perihelion_histograms_by_alignment(result, q_range, n_bins);
    let stats = population_statistics_series(result);

    let max_count = aligned
        .iter()
        .chain(anti_aligned.iter())
        .copied()
        .max()
        .unwrap_or(1)
        .max(1) as f64;

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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">{title}</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    let bin_width = plot_w / n_bins as f64;

    // Draw stacked bars: anti-aligned on bottom, aligned on top
    for (i, &center) in centers.iter().enumerate() {
        let x = margin + i as f64 * bin_width;
        let _ = center; // Used for positioning

        // Anti-aligned (blue)
        if anti_aligned[i] > 0 {
            let h = (anti_aligned[i] as f64 / max_count) * plot_h;
            let y = margin + plot_h - h;
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{:.1}" height="{h:.1}" fill="steelblue" stroke="black" stroke-width="0.3"/>"#,
                bin_width * 0.9,
            )
            .unwrap();
        }

        // Aligned (red) stacked on top
        if aligned[i] > 0 {
            let base_h = (anti_aligned[i] as f64 / max_count) * plot_h;
            let h = (aligned[i] as f64 / max_count) * plot_h;
            let y = margin + plot_h - base_h - h;
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{:.1}" height="{h:.1}" fill="indianred" stroke="black" stroke-width="0.3"/>"#,
                bin_width * 0.9,
            )
            .unwrap();
        }
    }

    // X-axis labels
    for q_val in [0, 50, 100, 150, 200, 250, 300, 350] {
        let x = margin + (q_val as f64 / q_range.1) * plot_w;
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{q_val}</text>"#,
            margin + plot_h + 15.0,
        )
        .unwrap();
    }
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Perihelion Distance q (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    // Y-axis label
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Count</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Legend
    write_legend(&mut svg, &stats, margin + plot_w - 140.0, margin + 15.0);

    writeln!(svg, "</svg>").unwrap();
    svg
}

fn write_legend(svg: &mut String, stats: &PopulationResult, x: f64, y: f64) {
    writeln!(
        svg,
        r#"<rect x="{x}" y="{y}" width="12" height="12" fill="indianred" stroke="black" stroke-width="0.3"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif">Aligned (n={})</text>"#,
        x + 16.0,
        y + 10.0,
        stats.aligned.count,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect x="{x}" y="{}" width="12" height="12" fill="steelblue" stroke="black" stroke-width="0.3"/>"#,
        y + 16.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif">Anti-aligned (n={})</text>"#,
        x + 16.0,
        y + 26.0,
        stats.anti_aligned.count,
    )
    .unwrap();
}

/// Generate an (a, q) scatter plot of the final snapshot colored by the
/// series-classified alignment (circulators in gray).
///
/// Useful for visualizing the two stable populations in the broad case.
pub fn aq_scatter_plot(result: &SimulationResult, title: &str, width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let a_range = (100.0, 600.0);
    let q_range = (0.0, 350.0);

    let alignments = classify_alignment_series(result);
    let final_snapshot = result.snapshots.last().expect("at least one snapshot");

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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">{title}</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Scatter points
    for (i, elem) in final_snapshot.elements.iter().enumerate() {
        let Some(elem) = elem else { continue };
        let q = elem.a * (1.0 - elem.e);
        let px = margin + (elem.a - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        let py = margin + plot_h - (q - q_range.0) / (q_range.1 - q_range.0) * plot_h;

        if px < margin || px > margin + plot_w || py < margin || py > margin + plot_h {
            continue;
        }

        let color = match alignments[i] {
            Some(Alignment::Aligned) => "indianred",
            Some(Alignment::AntiAligned) => "steelblue",
            Some(Alignment::Circulating) | None => "gray",
        };

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="3" fill="{color}" fill-opacity="0.6" stroke="black" stroke-width="0.3"/>"#,
        )
        .unwrap();
    }

    // Axis labels
    for a_val in [100, 200, 300, 400, 500, 600] {
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
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Semi-major Axis a (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    for q_val in [0, 50, 100, 150, 200, 250, 300, 350] {
        let y = margin + plot_h - (q_val as f64 - q_range.0) / (q_range.1 - q_range.0) * plot_h;
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="9" font-family="sans-serif">{q_val}</text>"#,
            margin - 5.0,
        )
        .unwrap();
    }
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Perihelion q (AU)</text>"#,
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
    use crate::simulation::KbSnapshot;
    use p9_core::types::OrbitalElements;
    use std::f64::consts::PI;

    /// Synthetic result: one aligned librator, one anti-aligned librator,
    /// one circulator (varpi_p9 = π).
    fn make_test_result() -> SimulationResult {
        let varpi_p9 = PI;
        let n_snaps = 12;
        let make_elem = |varpi: f64, a: f64, e: f64| {
            Some(OrbitalElements {
                a,
                e,
                i: 0.0,
                omega: varpi,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            })
        };
        let snapshots: Vec<KbSnapshot> = (0..n_snaps)
            .map(|k| {
                let phase = k as f64 * 0.7;
                let elements = vec![
                    // librates about ϖ = π (aligned with P9)
                    make_elem(PI + 0.4 * phase.sin(), 400.0, 0.7),
                    // librates about ϖ = 0 (anti-aligned)
                    make_elem(0.4 * phase.cos(), 300.0, 0.9),
                    // circulates
                    make_elem(k as f64 * 2.0 * PI / n_snaps as f64, 250.0, 0.8),
                ];
                KbSnapshot {
                    t: k as f64,
                    elements,
                    active_count: 3,
                    total_count: 3,
                }
            })
            .collect();
        SimulationResult {
            snapshots,
            config_summary: String::new(),
            varpi_p9,
        }
    }

    #[test]
    fn test_perihelion_distribution_plot() {
        let result = make_test_result();
        let svg = perihelion_distribution_plot(&result, "Test: Narrow Distribution", 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Perihelion Distance"));
        assert!(svg.contains("Aligned (n=1)"));
        assert!(svg.contains("Anti-aligned (n=1)"));
    }

    #[test]
    fn test_aq_scatter_plot() {
        let result = make_test_result();
        let svg = aq_scatter_plot(&result, "Test: a-q Scatter", 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Semi-major Axis"));
        assert!(svg.contains("circle"));
        // Circulator drawn in gray
        assert!(svg.contains("gray"));
    }
}
