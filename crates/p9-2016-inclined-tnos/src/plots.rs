//! SVG plot generation for inclined TNO figures.
//!
//! - Figure 1: density histograms in (a, i) and (q, i) space
//! - Scatter plots with known TNO overlays

use std::fmt::Write;

use p9_core::constants::*;

use crate::known_objects::KnownTno;
use crate::simulation::TnoSnapshot;

/// Generate a 2D density plot in (a, i) space.
pub fn density_plot_ai(
    snapshot: &TnoSnapshot,
    known_tnos: &[KnownTno],
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let a_range = (0.0, 600.0);
    let i_range = (0.0, 180.0);
    let n_a = 30;
    let n_i = 18;

    let (_a_bins, _i_bins, density) =
        crate::simulation::density_map_ai(snapshot, a_range, i_range, n_a, n_i);

    let max_count = density
        .iter()
        .flat_map(|row| row.iter())
        .cloned()
        .max()
        .unwrap_or(1)
        .max(1);

    let cell_w = plot_w / n_a as f64;
    let cell_h = plot_h / n_i as f64;

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
        r#"<text x="{}" y="18" text-anchor="middle" font-size="13" font-family="sans-serif">Particle Density (a, i)</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Density cells
    for (i_idx, row) in density.iter().enumerate() {
        for (a_idx, &count) in row.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let intensity = (count as f64 / max_count as f64).sqrt();
            let r = (255.0 * (1.0 - intensity)) as u8;
            let g = (255.0 * (1.0 - intensity)) as u8;
            let b = 255;

            let px = margin + a_idx as f64 * cell_w;
            let py = margin + (n_i - 1 - i_idx) as f64 * cell_h;

            writeln!(
                svg,
                r#"<rect x="{px:.1}" y="{py:.1}" width="{cell_w:.1}" height="{cell_h:.1}" fill="rgb({r},{g},{b})"/>"#,
            )
            .unwrap();
        }
    }

    // Known TNO markers
    for tno in known_tnos {
        let px = margin + (tno.elements.a - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        let i_deg = tno.elements.i * RAD2DEG;
        let py = margin + (1.0 - (i_deg - i_range.0) / (i_range.1 - i_range.0)) * plot_h;

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="5" fill="none" stroke="red" stroke-width="2"/>"#,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{:.1}" y="{:.1}" font-size="8" font-family="sans-serif" fill="red">{}</text>"#,
            px + 7.0,
            py + 3.0,
            tno.name,
        )
        .unwrap();
    }

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">a (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">i (°)</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    // 90° reference line
    let y90 = margin + (1.0 - 90.0 / 180.0) * plot_h;
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{y90:.1}" x2="{x2:.1}" y2="{y90:.1}" stroke="gray" stroke-dasharray="4,4"/>"#,
        x2 = margin + plot_w,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a scatter plot of surviving particles in (a, i) space.
pub fn scatter_plot_ai(
    snapshot: &TnoSnapshot,
    known_tnos: &[KnownTno],
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let a_range = (0.0, 600.0);
    let i_range = (0.0, 180.0);

    let mut svg = String::with_capacity(snapshot.elements.len() * 80 + 2000);
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

    for elem in &snapshot.elements {
        let i_deg = elem.i * RAD2DEG;
        if elem.a > a_range.1 || i_deg > i_range.1 {
            continue;
        }
        let px = margin + (elem.a - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        let py = margin + (1.0 - (i_deg - i_range.0) / (i_range.1 - i_range.0)) * plot_h;

        let color = if i_deg > 90.0 {
            "indianred"
        } else if i_deg > 50.0 {
            "orange"
        } else {
            "steelblue"
        };

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="2" fill="{color}" opacity="0.5"/>"#,
        )
        .unwrap();
    }

    // Known TNOs
    for tno in known_tnos {
        let px = margin + (tno.elements.a - a_range.0) / (a_range.1 - a_range.0) * plot_w;
        let i_deg = tno.elements.i * RAD2DEG;
        let py = margin + (1.0 - (i_deg - i_range.0) / (i_range.1 - i_range.0)) * plot_h;

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="6" fill="none" stroke="black" stroke-width="2"/>"#,
        )
        .unwrap();
    }

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">a (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">i (°)</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_plot() {
        let snapshot = TnoSnapshot {
            t: 0.0,
            elements: vec![
                p9_core::types::OrbitalElements {
                    a: 50.0,
                    e: 0.5,
                    i: 100.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                p9_core::types::OrbitalElements {
                    a: 300.0,
                    e: 0.7,
                    i: 30.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 2,
            total_count: 2,
        };

        let known = crate::known_objects::paper_tnos();
        let svg = density_plot_ai(&snapshot, &known, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Drac"));
    }

    #[test]
    fn test_scatter_plot() {
        let snapshot = TnoSnapshot {
            t: 0.0,
            elements: vec![p9_core::types::OrbitalElements {
                a: 50.0,
                e: 0.5,
                i: 100.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            }],
            active_count: 1,
            total_count: 1,
        };

        let known = crate::known_objects::paper_tnos();
        let svg = scatter_plot_ai(&snapshot, &known, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("circle"));
    }
}
