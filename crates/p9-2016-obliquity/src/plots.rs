//! SVG plot generation for obliquity figures.
//!
//! - Figure 2: contours of required i₉ in (a₉, e₉) plane
//! - Time evolution of solar obliquity over 4.5 Gyr

use std::fmt::Write;

use p9_core::constants::*;

use crate::parameter_survey::SurveyResult;
use crate::secular_hamiltonian::ObliquitySnapshot;

/// Generate a contour-like plot of required i₉ in (a₉, e₉) space.
pub fn required_inclination_plot(
    results: &[SurveyResult],
    mass_earth: f64,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let filtered: Vec<&SurveyResult> = results
        .iter()
        .filter(|r| (r.mass_earth - mass_earth).abs() < 0.01)
        .collect();

    let mut svg = String::with_capacity(filtered.len() * 100 + 1000);
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
        r#"<text x="{}" y="18" text-anchor="middle" font-size="13" font-family="sans-serif">Required i₉ for 6° obliquity (m = {} M⊕)</text>"#,
        width as f64 / 2.0,
        mass_earth,
    )
    .unwrap();

    let a_min = 350.0_f64;
    let a_max = 950.0_f64;
    let e_min = 0.2_f64;
    let e_max = 0.95_f64;

    let cell_w = plot_w / 7.0;
    let cell_h = plot_h / 7.0;

    for result in &filtered {
        if let Some(i9) = result.required_i9 {
            let i9_deg = i9 * RAD2DEG;
            let px = margin + (result.a9 - a_min) / (a_max - a_min) * plot_w;
            let py = margin + (1.0 - (result.e9 - e_min) / (e_max - e_min)) * plot_h;

            // Color by required inclination: blue (low) to red (high)
            let t = ((i9_deg - 10.0) / 30.0).clamp(0.0, 1.0);
            let r = (255.0 * t) as u8;
            let b = (255.0 * (1.0 - t)) as u8;

            writeln!(
                svg,
                r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="rgb({r},50,{b})" opacity="0.7"/>"#,
                px - cell_w / 2.0,
                py - cell_h / 2.0,
                cell_w,
                cell_h,
            )
            .unwrap();

            // Label with required i₉
            writeln!(
                svg,
                r#"<text x="{px:.1}" y="{:.1}" text-anchor="middle" font-size="9" font-family="sans-serif" fill="white">{:.0}°</text>"#,
                py + 3.0,
                i9_deg,
            )
            .unwrap();
        }
    }

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">a₉ (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">e₉</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a time evolution plot of solar obliquity.
pub fn obliquity_evolution_plot(
    snapshots: &[ObliquitySnapshot],
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let t_max = snapshots.last().map_or(1.0, |s| s.t);
    let obl_max = snapshots
        .iter()
        .map(|s| s.obliquity * RAD2DEG)
        .fold(0.0_f64, f64::max)
        .max(10.0);

    let mut svg = String::with_capacity(snapshots.len() * 60 + 1000);
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
        r#"<text x="{}" y="18" text-anchor="middle" font-size="13" font-family="sans-serif">Solar Obliquity Evolution</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot obliquity curve
    let mut path = String::from("M");
    for (i, snap) in snapshots.iter().enumerate() {
        let px = margin + (snap.t / t_max) * plot_w;
        let py = margin + (1.0 - snap.obliquity * RAD2DEG / obl_max) * plot_h;
        if i == 0 {
            write!(path, "{px:.1},{py:.1}").unwrap();
        } else {
            write!(path, " L{px:.1},{py:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<path d="{path}" fill="none" stroke="steelblue" stroke-width="2"/>"#,
    )
    .unwrap();

    // 6° reference line
    let ref_y = margin + (1.0 - 6.0 / obl_max) * plot_h;
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{ref_y:.1}" x2="{x2:.1}" y2="{ref_y:.1}" stroke="red" stroke-dasharray="4,4"/>"#,
        x2 = margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{x:.1}" y="{y:.1}" font-size="9" font-family="sans-serif" fill="red">6° observed</text>"#,
        x = margin + plot_w + 3.0,
        y = ref_y + 3.0,
    )
    .unwrap();

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Time (Gyr)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">Obliquity (°)</text>"#,
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
    fn test_required_inclination_plot() {
        let results = vec![SurveyResult {
            mass_earth: 10.0,
            a9: 700.0,
            e9: 0.6,
            required_i9: Some(25.0 * DEG2RAD),
            final_obliquity_deg: 6.0,
        }];

        let svg = required_inclination_plot(&results, 10.0, 500, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_obliquity_evolution_plot() {
        let snapshots = vec![
            ObliquitySnapshot {
                t: 0.0,
                obliquity: 0.001,
                mutual_inclination: 0.35,
                omega_big_sun: 0.0,
                i_9: 0.35,
                omega_big_9: 3.14,
                delta_omega_big: -3.14,
            },
            ObliquitySnapshot {
                t: 4.5 * GYR_DAYS,
                obliquity: 6.0 * DEG2RAD,
                mutual_inclination: 0.34,
                omega_big_sun: 1.2,
                i_9: 0.34,
                omega_big_9: 2.5,
                delta_omega_big: -1.3,
            },
        ];

        let svg = obliquity_evolution_plot(&snapshots, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("path"));
        assert!(svg.contains("6°"));
    }
}
