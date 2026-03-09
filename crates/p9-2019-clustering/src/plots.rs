//! SVG plot generation for clustering analysis figures.
//!
//! Reproduces Figures 1-2 from Brown & Batygin (2019):
//! - (x, y) Poincaré plane with observed and Monte Carlo points
//! - (p, q) pole position plane

use std::fmt::Write;

use crate::clustering_analysis::ClusteringResult;
use crate::poincare_variables::PoincareState;

/// Generate a Poincaré (x, y) clustering plot.
///
/// Shows individual KBO positions and the mean clustering point,
/// with p-value annotation.
pub fn perihelion_clustering_plot(
    states: &[PoincareState],
    result: &ClusteringResult,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    // Data range centered on origin
    let range = 1.2_f64;
    let cx = margin + plot_w / 2.0;
    let cy = margin + plot_h / 2.0;
    let scale = plot_w.min(plot_h) / (2.0 * range);

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Perihelion Clustering (p = {:.1}%)</text>"#,
        width as f64 / 2.0,
        result.p_perihelion * 100.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Axes through origin
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{cy:.1}" x2="{}" y2="{cy:.1}" stroke="gray" stroke-width="0.5" stroke-dasharray="4,4"/>"#,
        margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<line x1="{cx:.1}" y1="{margin}" x2="{cx:.1}" y2="{}" stroke="gray" stroke-width="0.5" stroke-dasharray="4,4"/>"#,
        margin + plot_h,
    )
    .unwrap();

    // KBO points
    for state in states {
        let px = cx + state.x * scale;
        let py = cy - state.y * scale;
        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="4" fill="steelblue" stroke="black" stroke-width="0.5"/>"#,
        )
        .unwrap();
    }

    // Mean point
    let mean = crate::poincare_variables::mean_state(states);
    let mx = cx + mean.x * scale;
    let my = cy - mean.y * scale;
    writeln!(
        svg,
        r#"<circle cx="{mx:.1}" cy="{my:.1}" r="6" fill="red" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    // Origin
    writeln!(
        svg,
        r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="3" fill="none" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">x = sqrt(2Γ) cos(ϖ)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 30.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">y = sqrt(2Γ) sin(ϖ)</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a pole position (p, q) clustering plot.
pub fn pole_clustering_plot(
    states: &[PoincareState],
    result: &ClusteringResult,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let range = 0.6_f64;
    let cx = margin + plot_w / 2.0;
    let cy = margin + plot_h / 2.0;
    let scale = plot_w.min(plot_h) / (2.0 * range);

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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Pole Clustering (p = {:.1}%)</text>"#,
        width as f64 / 2.0,
        result.p_pole * 100.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Axes
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{cy:.1}" x2="{}" y2="{cy:.1}" stroke="gray" stroke-width="0.5" stroke-dasharray="4,4"/>"#,
        margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<line x1="{cx:.1}" y1="{margin}" x2="{cx:.1}" y2="{}" stroke="gray" stroke-width="0.5" stroke-dasharray="4,4"/>"#,
        margin + plot_h,
    )
    .unwrap();

    // KBO pole positions
    for state in states {
        let px = cx + state.p * scale;
        let py = cy - state.q_var * scale;
        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="4" fill="steelblue" stroke="black" stroke-width="0.5"/>"#,
        )
        .unwrap();
    }

    // Mean point
    let mean = crate::poincare_variables::mean_state(states);
    let mx = cx + mean.p * scale;
    let my = cy - mean.q_var * scale;
    writeln!(
        svg,
        r#"<circle cx="{mx:.1}" cy="{my:.1}" r="6" fill="red" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">p = sqrt(2Z) cos(Ω)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 30.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">q = sqrt(2Z) sin(Ω)</text>"#,
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
    use crate::clustering_analysis;
    use crate::kbo_sample;

    #[test]
    fn test_perihelion_plot() {
        let kbos = kbo_sample::paper_sample_a230();
        let states = clustering_analysis::compute_poincare_states(&kbos);
        let result = ClusteringResult {
            observed_perihelion: 0.45,
            observed_pole: 0.13,
            observed_combined: 0.47,
            mean_varpi: 1.27,
            mean_omega: 1.43,
            p_perihelion: 0.04,
            p_pole: 0.035,
            p_combined: 0.002,
            n_iterations: 10000,
        };

        let svg = perihelion_clustering_plot(&states, &result, 500, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Perihelion Clustering"));
    }

    #[test]
    fn test_pole_plot() {
        let kbos = kbo_sample::paper_sample_a230();
        let states = clustering_analysis::compute_poincare_states(&kbos);
        let result = ClusteringResult {
            observed_perihelion: 0.45,
            observed_pole: 0.13,
            observed_combined: 0.47,
            mean_varpi: 1.27,
            mean_omega: 1.43,
            p_perihelion: 0.04,
            p_pole: 0.035,
            p_combined: 0.002,
            n_iterations: 10000,
        };

        let svg = pole_clustering_plot(&states, &result, 500, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Pole Clustering"));
    }
}
