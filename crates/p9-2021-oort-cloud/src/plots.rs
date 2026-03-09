//! SVG plot generation for the Oort cloud injection paper figures.
//!
//! Produces:
//! - Confinement comparison bar chart (f_varpi for IOC vs scattered disk)
//! - Semi-major axis distribution stacked histogram

use std::fmt::Write;

use crate::injection_simulation::InjectionResult;
use crate::population_comparison::{semi_major_axis_distribution, PopulationComparison};

/// Generate a bar chart comparing f_varpi confinement between populations.
///
/// The IOC population shows weaker confinement (~67%) compared to the
/// scattered disk (~88%), visualized as side-by-side bars.
pub fn confinement_comparison_plot(
    comparison: &PopulationComparison,
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="13" font-family="sans-serif">Longitude of Perihelion Confinement (f_varpi)</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot area border
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    // Y-axis scale: 0 to 1
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = margin + plot_h * (1.0 - tick);
        writeln!(
            svg,
            "<line x1=\"{margin}\" y1=\"{y:.1}\" x2=\"{}\" y2=\"{y:.1}\" stroke=\"#ddd\" stroke-width=\"1\"/>",
            margin + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{:.1}" text-anchor="end" font-size="9" font-family="sans-serif">{:.0}%</text>"#,
            margin - 5.0,
            y + 3.0,
            tick * 100.0,
        )
        .unwrap();
    }

    // Bar width and positions
    let bar_w = plot_w * 0.25;
    let gap = plot_w * 0.1;
    let x_ioc = margin + plot_w / 2.0 - bar_w - gap / 2.0;
    let x_sd = margin + plot_w / 2.0 + gap / 2.0;

    // IOC bar
    let h_ioc = plot_h * comparison.f_varpi_ioc;
    let y_ioc = margin + plot_h - h_ioc;
    writeln!(
        svg,
        r#"<rect x="{x_ioc:.1}" y="{y_ioc:.1}" width="{bar_w:.1}" height="{h_ioc:.1}" fill="steelblue" opacity="0.8"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="10" font-family="sans-serif" font-weight="bold">{:.0}%</text>"#,
        x_ioc + bar_w / 2.0,
        y_ioc - 5.0,
        comparison.f_varpi_ioc * 100.0,
    )
    .unwrap();

    // Scattered disk bar
    let h_sd = plot_h * comparison.f_varpi_scattered;
    let y_sd = margin + plot_h - h_sd;
    writeln!(
        svg,
        r#"<rect x="{x_sd:.1}" y="{y_sd:.1}" width="{bar_w:.1}" height="{h_sd:.1}" fill="indianred" opacity="0.8"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="10" font-family="sans-serif" font-weight="bold">{:.0}%</text>"#,
        x_sd + bar_w / 2.0,
        y_sd - 5.0,
        comparison.f_varpi_scattered * 100.0,
    )
    .unwrap();

    // Labels
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{}" text-anchor="middle" font-size="10" font-family="sans-serif">IOC</text>"#,
        x_ioc + bar_w / 2.0,
        margin + plot_h + 15.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{}" text-anchor="middle" font-size="10" font-family="sans-serif">Scattered Disk</text>"#,
        x_sd + bar_w / 2.0,
        margin + plot_h + 15.0,
    )
    .unwrap();

    // Y-axis label
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">f_varpi</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a stacked histogram of the semi-major axis distribution.
///
/// Shows IOC-injected objects (blue) stacked with scattered disk objects (red),
/// highlighting the IOC enhancement at a > 2000 AU.
pub fn sma_distribution_plot(result: &InjectionResult, width: u32, height: u32) -> String {
    let bins = semi_major_axis_distribution(result, 20);

    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let max_count = bins
        .iter()
        .map(|b| b.count_ioc + b.count_scattered)
        .max()
        .unwrap_or(1)
        .max(1);

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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="13" font-family="sans-serif">Semi-Major Axis Distribution: IOC vs Scattered Disk</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot area border
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    let n_bins = bins.len();
    let bar_w = plot_w / n_bins as f64;
    let a_min = bins.first().map(|b| b.a_low).unwrap_or(250.0);
    let a_max = bins.last().map(|b| b.a_high).unwrap_or(20000.0);

    // Draw stacked bars
    for (i, bin) in bins.iter().enumerate() {
        let x = margin + i as f64 * bar_w;

        // Scattered disk (bottom, red)
        let h_sd = plot_h * bin.count_scattered as f64 / max_count as f64;
        let y_sd = margin + plot_h - h_sd;
        if bin.count_scattered > 0 {
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y_sd:.1}" width="{bar_w:.1}" height="{h_sd:.1}" fill="indianred" opacity="0.7"/>"#,
            )
            .unwrap();
        }

        // IOC (top, blue)
        let h_ioc = plot_h * bin.count_ioc as f64 / max_count as f64;
        let y_ioc = y_sd - h_ioc;
        if bin.count_ioc > 0 {
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y_ioc:.1}" width="{bar_w:.1}" height="{h_ioc:.1}" fill="steelblue" opacity="0.7"/>"#,
            )
            .unwrap();
        }
    }

    // X-axis tick labels (select a few representative values)
    let tick_values = [500.0, 2000.0, 5000.0, 10000.0, 15000.0, 20000.0];
    for &a_val in &tick_values {
        if a_val >= a_min && a_val <= a_max {
            let x = margin + (a_val - a_min) / (a_max - a_min) * plot_w;
            writeln!(
                svg,
                r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="8" font-family="sans-serif">{}</text>"#,
                margin + plot_h + 15.0,
                a_val as u64,
            )
            .unwrap();
        }
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Semi-Major Axis a (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
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
    let lx = margin + plot_w - 120.0;
    let ly = margin + 15.0;
    writeln!(
        svg,
        r#"<rect x="{lx:.1}" y="{ly:.1}" width="12" height="12" fill="steelblue" opacity="0.7"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="9" font-family="sans-serif">IOC Injected</text>"#,
        lx + 16.0,
        ly + 10.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect x="{lx:.1}" y="{:.1}" width="12" height="12" fill="indianred" opacity="0.7"/>"#,
        ly + 18.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="9" font-family="sans-serif">Scattered Disk</text>"#,
        lx + 16.0,
        ly + 28.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_comparison() -> PopulationComparison {
        PopulationComparison {
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            n_ioc_injected: 50,
            n_ioc_total: 300,
            injection_efficiency: 50.0 / 300.0,
        }
    }

    fn sample_result() -> InjectionResult {
        InjectionResult {
            injection_fraction: 0.1,
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            n_injected: 5,
            n_total: 50,
            injected_sma: vec![500.0, 3000.0, 5000.0, 8000.0, 15000.0],
            injected_dvarpi: vec![0.1, -0.5, 0.3, -0.2, 0.7],
        }
    }

    #[test]
    fn test_confinement_comparison_plot_valid_svg() {
        let svg = confinement_comparison_plot(&sample_comparison(), 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("IOC"));
        assert!(svg.contains("Scattered Disk"));
        assert!(svg.contains("f_varpi"));
    }

    #[test]
    fn test_sma_distribution_plot_valid_svg() {
        let svg = sma_distribution_plot(&sample_result(), 700, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Semi-Major Axis"));
        assert!(svg.contains("IOC Injected"));
        assert!(svg.contains("steelblue"));
        assert!(svg.contains("indianred"));
    }
}
