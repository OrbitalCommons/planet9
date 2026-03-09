//! SVG plot generation for paper figures.
//!
//! Generates SVG visualizations reproducing the key figures from
//! Batygin & Brown (2016):
//! - Figure 1: KBO orbital element clustering
//! - Figure 3/4: Phase-space portraits (analytical + N-body)
//! - Figure 5: Scattered disk Δϖ footprint
//! - Figure 8: 3D scattered disk clustering

use std::fmt::Write;

use crate::kbo_elements::KboRecord;
use crate::phase_portrait::TrajectoryPoint;
use crate::scattered_disk_sim::DiskSnapshot;
use p9_core::constants::*;

/// Generate Figure 1: KBO orbital angles vs semi-major axis.
///
/// Three panels: ϖ vs a, Ω vs a, ω vs a
pub fn figure1_kbo_clustering(kbos: &[KboRecord], width: u32, height: u32) -> String {
    let panel_h = height / 3;
    let margin = 50.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = panel_h as f64 - 2.0 * margin / 3.0;

    let mut svg = String::with_capacity(5000);
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

    let a_min = 100.0_f64;
    let a_max = 600.0_f64;
    let a_range = a_max - a_min;

    // Panel labels and data extractors
    let panels: Vec<(&str, Box<dyn Fn(&KboRecord) -> f64>)> = vec![
        (
            "ϖ (deg)",
            Box::new(|k: &KboRecord| (k.elements.omega + k.elements.omega_big) * RAD2DEG),
        ),
        (
            "Ω (deg)",
            Box::new(|k: &KboRecord| k.elements.omega_big * RAD2DEG),
        ),
        (
            "ω (deg)",
            Box::new(|k: &KboRecord| k.elements.omega * RAD2DEG),
        ),
    ];

    for (panel_idx, (label, extractor)) in panels.iter().enumerate() {
        let y_offset = panel_idx as f64 * panel_h as f64;

        // Panel border
        writeln!(
            svg,
            r#"<rect x="{margin}" y="{y:.1}" width="{plot_w}" height="{plot_h:.1}" fill="none" stroke="black"/>"#,
            y = y_offset + margin / 3.0,
        )
        .unwrap();

        // Y-axis label
        writeln!(
            svg,
            r#"<text x="12" y="{y:.1}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90,12,{y:.1})">{label}</text>"#,
            y = y_offset + margin / 3.0 + plot_h / 2.0,
        )
        .unwrap();

        // Plot KBO points
        for kbo in kbos {
            let angle = extractor(kbo);
            let px = margin + (kbo.elements.a - a_min) / a_range * plot_w;
            let py = y_offset + margin / 3.0 + (1.0 - angle / 360.0) * plot_h;

            writeln!(
                svg,
                r#"<circle cx="{px:.1}" cy="{py:.1}" r="5" fill="steelblue" stroke="black" stroke-width="0.5"/>"#,
            )
            .unwrap();
        }
    }

    // X-axis label
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">a (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 5.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a phase-space portrait figure (e vs Δϖ).
///
/// Can be used for both analytical (contour) and N-body (scatter) portraits.
pub fn phase_portrait_figure(
    points: &[TrajectoryPoint],
    title: &str,
    width: u32,
    height: u32,
) -> String {
    let margin = 50.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(points.len() * 80 + 1000);
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
        r#"<text x="{}" y="15" text-anchor="middle" font-size="13" font-family="sans-serif">{title}</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot points
    for p in points {
        let px = margin + (p.delta_varpi + std::f64::consts::PI) / TWO_PI * plot_w;
        let py = margin + (1.0 - p.e / 0.95) * plot_h;

        // Color by alignment: anti-aligned (|Δϖ|>π/2) = blue, aligned = red
        let color = if p.delta_varpi.abs() > std::f64::consts::PI / 2.0 {
            "steelblue"
        } else {
            "indianred"
        };

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="1" fill="{color}" opacity="0.4"/>"#,
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
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Δϖ (rad)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">e</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a scattered disk clustering figure (ϖ, Ω, ω vs a).
pub fn scattered_disk_clustering(
    snapshot: &DiskSnapshot,
    varpi_p9: f64,
    width: u32,
    height: u32,
) -> String {
    let margin = 50.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(snapshot.elements.len() * 100 + 1000);
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

    let a_min = 50.0_f64;
    let a_max = 600.0_f64;
    let a_range = a_max - a_min;

    for elem in &snapshot.elements {
        if elem.a < a_min || elem.a > a_max {
            continue;
        }

        let varpi = elem.omega + elem.omega_big;
        let mut dv = varpi - varpi_p9;
        while dv > std::f64::consts::PI {
            dv -= TWO_PI;
        }
        while dv < -std::f64::consts::PI {
            dv += TWO_PI;
        }

        let px = margin + (elem.a - a_min) / a_range * plot_w;
        let py = margin + (1.0 - (dv + std::f64::consts::PI) / TWO_PI) * plot_h;

        let color = if dv.abs() > std::f64::consts::PI / 2.0 {
            "steelblue"
        } else {
            "indianred"
        };

        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="3" fill="{color}" opacity="0.7"/>"#,
        )
        .unwrap();
    }

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
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">Δϖ (rad)</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kbo_elements::stable_kbos;
    use p9_core::types::OrbitalElements;

    #[test]
    fn test_figure1_generation() {
        let kbos = stable_kbos();
        let svg = figure1_kbo_clustering(&kbos, 600, 600);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("circle"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_phase_portrait_figure() {
        let points = vec![
            TrajectoryPoint {
                t: 0.0,
                e: 0.3,
                delta_varpi: 0.5,
                a: 300.0,
            },
            TrajectoryPoint {
                t: 0.0,
                e: 0.7,
                delta_varpi: -2.0,
                a: 300.0,
            },
        ];
        let svg = phase_portrait_figure(&points, "a = 300 AU", 400, 300);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("300 AU"));
    }

    #[test]
    fn test_scattered_disk_figure() {
        let snapshot = DiskSnapshot {
            t: 0.0,
            elements: vec![
                OrbitalElements {
                    a: 300.0,
                    e: 0.7,
                    i: 0.0,
                    omega: 5.5,
                    omega_big: 0.5,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 400.0,
                    e: 0.8,
                    i: 0.1,
                    omega: 5.0,
                    omega_big: 1.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 2,
            total_count: 2,
        };
        let svg = scattered_disk_clustering(&snapshot, 2.5, 500, 400);
        assert!(svg.contains("<svg"));
    }
}
