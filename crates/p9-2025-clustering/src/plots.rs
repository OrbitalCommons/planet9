//! SVG plot generation for clustering analysis figures.

use std::fmt::Write;

use crate::clustering::{von_mises_pdf, VonMisesParams};
use crate::orbital_poles::OrbitalPole;

/// Plot longitude of perihelion distribution with von Mises fit overlay.
pub fn varpi_distribution_plot(
    angles: &[f64],
    fit: &VonMisesParams,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(6000);
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Longitude of Perihelion Clustering</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Histogram of angles
    let n_bins = 18; // 20-degree bins
    let mut bins = vec![0usize; n_bins];
    let bin_width_rad = std::f64::consts::TAU / n_bins as f64;
    for &a in angles {
        let a_norm = ((a % std::f64::consts::TAU) + std::f64::consts::TAU) % std::f64::consts::TAU;
        let idx = (a_norm / bin_width_rad) as usize;
        let idx = idx.min(n_bins - 1);
        bins[idx] += 1;
    }

    let max_count = *bins.iter().max().unwrap_or(&1) as f64;
    let bar_w = plot_w / n_bins as f64;

    for (i, &count) in bins.iter().enumerate() {
        let h = (count as f64 / max_count) * plot_h * 0.8;
        let x = margin + i as f64 * bar_w;
        let y = margin + plot_h - h;
        writeln!(
            svg,
            r#"<rect x="{x:.1}" y="{y:.1}" width="{:.1}" height="{h:.1}" fill="steelblue" opacity="0.6"/>"#,
            bar_w * 0.9,
        )
        .unwrap();
    }

    // Von Mises fit curve
    let n_curve = 200;
    let mut curve_points = Vec::with_capacity(n_curve);
    let mut pdf_max = 0.0_f64;
    for ic in 0..n_curve {
        let theta = ic as f64 * std::f64::consts::TAU / n_curve as f64;
        let pdf = von_mises_pdf(theta, fit.mu, fit.kappa);
        pdf_max = pdf_max.max(pdf);
    }

    for ic in 0..n_curve {
        let theta = ic as f64 * std::f64::consts::TAU / n_curve as f64;
        let pdf = von_mises_pdf(theta, fit.mu, fit.kappa);
        let x = margin + (theta / std::f64::consts::TAU) * plot_w;
        let y = margin + plot_h - (pdf / pdf_max.max(1e-10)) * plot_h * 0.8;
        curve_points.push(format!("{x:.1},{y:.1}"));
    }

    if !curve_points.is_empty() {
        writeln!(
            svg,
            r#"<polyline points="{}" fill="none" stroke="crimson" stroke-width="2"/>"#,
            curve_points.join(" "),
        )
        .unwrap();
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Longitude of perihelion (deg)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Plot orbital pole positions on the (ix, iy) plane.
pub fn orbital_pole_plot(poles: &[OrbitalPole], width: u32, height: u32) -> String {
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

    writeln!(
        svg,
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Orbital Pole Distribution</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    if !poles.is_empty() {
        let max_r = poles
            .iter()
            .map(|p| p.inclination())
            .fold(0.0_f64, f64::max)
            .max(0.01);

        let cx = margin + plot_w / 2.0;
        let cy = margin + plot_h / 2.0;
        let scale = (plot_w.min(plot_h) / 2.0) / max_r;

        for pole in poles {
            let px = cx + pole.x * scale;
            let py = cy - pole.y * scale;
            writeln!(
                svg,
                r#"<circle cx="{px:.1}" cy="{py:.1}" r="4" fill="steelblue" opacity="0.7"/>"#,
            )
            .unwrap();
        }
    }

    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">i cos(Omega)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varpi_plot_valid_svg() {
        let angles = vec![0.5, 0.8, 0.9, 1.0, 1.1, 0.7];
        let fit = VonMisesParams {
            mu: 0.87,
            kappa: 5.0,
        };
        let svg = varpi_distribution_plot(&angles, &fit, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_pole_plot_valid_svg() {
        let poles = vec![
            OrbitalPole {
                name: "a",
                x: 0.1,
                y: 0.2,
            },
            OrbitalPole {
                name: "b",
                x: -0.1,
                y: 0.15,
            },
        ];
        let svg = orbital_pole_plot(&poles, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_varpi_plot_empty() {
        let fit = VonMisesParams {
            mu: 0.0,
            kappa: 1.0,
        };
        let svg = varpi_distribution_plot(&[], &fit, 600, 400);
        assert!(svg.contains("<svg"));
    }
}
