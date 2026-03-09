//! SVG plot generation for bias and clustering figures.

use std::fmt::Write;

use crate::clustering_test::ClusteringResult;
use crate::kbo_sample::DistantKbo;

/// Generate a plot showing KBO ϖ distribution with clustering statistic.
pub fn varpi_distribution_plot(
    kbos: &[DistantKbo],
    result: &ClusteringResult,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(2000);
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
        r#"<text x="{}" y="18" text-anchor="middle" font-size="13" font-family="sans-serif">ϖ Distribution (p = {:.3}%)</text>"#,
        width as f64 / 2.0,
        result.p_varpi * 100.0,
    )
    .unwrap();

    // Plot KBO ϖ values as points on a unit circle projection
    let cx = margin + plot_w / 2.0;
    let cy = margin + plot_h / 2.0;
    let r = plot_w.min(plot_h) / 2.5;

    // Reference circle
    writeln!(
        svg,
        r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="none" stroke="gray"/>"#,
    )
    .unwrap();

    // KBO points
    for kbo in kbos {
        let varpi = kbo.elements.omega + kbo.elements.omega_big;
        let px = cx + r * varpi.cos();
        let py = cy - r * varpi.sin();
        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="5" fill="steelblue" stroke="black" stroke-width="0.5"/>"#,
        )
        .unwrap();
    }

    // Mean direction
    let varpis: Vec<f64> = kbos
        .iter()
        .map(|k| k.elements.omega + k.elements.omega_big)
        .collect();
    let sin_sum: f64 = varpis.iter().map(|v| v.sin()).sum();
    let cos_sum: f64 = varpis.iter().map(|v| v.cos()).sum();
    let mean_varpi = sin_sum.atan2(cos_sum);
    let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / kbos.len() as f64;

    let mx = cx + r * r_bar * mean_varpi.cos();
    let my = cy - r * r_bar * mean_varpi.sin();
    writeln!(
        svg,
        r#"<line x1="{cx:.1}" y1="{cy:.1}" x2="{mx:.1}" y2="{my:.1}" stroke="red" stroke-width="2"/>"#,
    )
    .unwrap();

    // Axes labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="10" font-family="sans-serif">0°</text>"#,
        cx + r + 15.0,
        cy + 4.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="10" font-family="sans-serif">180°</text>"#,
        cx - r - 15.0,
        cy + 4.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kbo_sample;

    #[test]
    fn test_varpi_distribution_plot() {
        let kbos = kbo_sample::paper_sample_a230();
        let result = ClusteringResult {
            observed_z_varpi: 5.0,
            observed_z_omega: 3.0,
            p_varpi: 0.012,
            p_omega: 0.05,
            p_combined: 0.00025,
            n_iterations: 10000,
        };

        let svg = varpi_distribution_plot(&kbos, &result, 500, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("circle"));
    }
}
