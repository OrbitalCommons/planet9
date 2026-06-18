//! Reproduction of Linder & Mordasini (2016), "Evolution and magnitudes of
//! candidate Planet Nine" (arXiv:1602.07465).
//!
//! The paper couples a planetary thermal-evolution model (a cooling, several-
//! Earth-mass ice/gas planet) to a grey atmosphere to predict Planet Nine's
//! brightness across photometric bands as a function of mass and heliocentric
//! distance. The headline result is twofold:
//!
//! 1. In **reflected sunlight** (optical/near-IR: V, R, I) the planet is faint
//!    — a 10 M⊕ body at ~700 AU has V ≈ 21.7 — and gets rapidly fainter with
//!    decreasing mass and increasing distance, because reflected flux scales as
//!    R²/(r²Δ²) (a 1/r⁴-like falloff at large r).
//!
//! 2. In **thermal emission** the planet radiates its internal heat at an
//!    effective temperature of ~40–50 K. On the steep Wien tail this makes the
//!    near-IR (WISE W1/W2 ≈ 3–5 µm) hostile but the **mid/far-IR longward of
//!    ~13 µm** (the N, Q and WISE W4 bands) far more favorable: the paper finds
//!    Planet Nine is brightest and most detectable in the thermal mid-IR.
//!
//! This crate computes both channels from real physics:
//!
//! - The reflected-light V magnitude reuses
//!   [`p9_core::analysis::photometry::planet_apparent_magnitude`] (the shared
//!   Neptune-anchored mass-radius relation + the reflected-sunlight apparent-
//!   magnitude law), so it cannot drift from the rest of the workspace.
//!
//! - The thermal magnitudes come from a grey-body Planck emitter at the
//!   evolution-model effective temperature, integrated against each band's
//!   zero point ([`thermal`]).
//!
//! Published predictions are kept as labelled reference constants
//! ([`reference`]) and the residuals are documented honestly: the reflected
//! V-band model and the single-temperature grey-body IR model are both
//! simplifications of the paper's coupled evolution+atmosphere grid.

pub mod reference;
pub mod thermal;

use p9_core::analysis::photometry::planet_apparent_magnitude;
use p9_core::units::{au, earth_masses, Length, Mass};

/// Neptune-like geometric albedo used for the reflected-light channel, matching
/// the workspace default (`p9_core::analysis::photometry::ALBEDO_NEPTUNE`).
pub const P9_ALBEDO: f64 = p9_core::analysis::photometry::ALBEDO_NEPTUNE;

/// Predicted reflected-light apparent V magnitude of a Planet Nine of
/// `mass_earth` Earth masses at heliocentric distance `distance_au`, observed
/// at opposition. Thin wrapper over the shared workspace photometry.
pub fn reflected_v_magnitude(mass_earth: f64, distance_au: f64) -> f64 {
    planet_apparent_magnitude(mass_earth, P9_ALBEDO, distance_au)
}

/// One row of the magnitude-vs-mass table at a fixed distance: the reflected V
/// magnitude and the thermal magnitude in the most-favorable mid-IR band
/// (WISE W4, 22 µm), for a planet of `mass_earth` at `distance_au`.
#[derive(Debug, Clone, Copy)]
pub struct MagnitudeRow {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
    /// Effective temperature of the cooling planet (K), from the evolution
    /// model floor (see [`thermal::effective_temperature`]).
    pub t_eff_k: f64,
    /// Reflected-light apparent V magnitude.
    pub v_mag: f64,
    /// Thermal apparent magnitude in the most-favorable mid-IR band (WISE W4).
    pub mid_ir_mag: f64,
}

impl MagnitudeRow {
    // ---- typed (uom) boundary: dimension-checked views of the f64 fields ----

    /// Planet mass as a [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }
    /// Heliocentric distance as a [`Length`].
    pub fn distance(&self) -> Length {
        au(self.distance_au)
    }
}

/// Build a [`MagnitudeRow`] for the given mass and distance.
pub fn magnitude_row(mass_earth: f64, distance_au: f64) -> MagnitudeRow {
    let planet = thermal::P9 {
        mass_earth,
        distance_au,
    };
    MagnitudeRow {
        mass_earth,
        distance_au,
        t_eff_k: planet.effective_temperature(),
        v_mag: reflected_v_magnitude(mass_earth, distance_au),
        mid_ir_mag: planet.thermal_magnitude(thermal::Band::W4),
    }
}

/// Tabulate [`MagnitudeRow`]s for a set of masses at a fixed distance — the
/// paper's "magnitudes vs mass" presentation.
pub fn magnitude_table(masses_earth: &[f64], distance_au: f64) -> Vec<MagnitudeRow> {
    masses_earth
        .iter()
        .map(|&m| magnitude_row(m, distance_au))
        .collect()
}

/// Whether a body of apparent magnitude `mag` in band `band` is brighter than
/// (i.e. detectable by) a survey with limiting magnitude `limiting_mag`.
pub fn is_detectable(mag: f64, limiting_mag: f64) -> bool {
    mag <= limiting_mag
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::surveys::limiting_magnitude;
    use p9_core::units::EARTH_MASS_KG;
    use uom::si::length::astronomical_unit;
    use uom::si::mass::kilogram;

    #[test]
    fn typed_row_accessors_match_f64_fields() {
        let row = magnitude_row(10.0, 700.0);
        assert_relative_eq!(
            row.mass().get::<kilogram>(),
            row.mass_earth * EARTH_MASS_KG,
            epsilon = 1.0
        );
        assert_relative_eq!(row.distance().get::<astronomical_unit>(), row.distance_au);
    }

    #[test]
    fn ten_earth_mass_at_700au_matches_published_v() {
        // Headline of the paper: a 10 M⊕ Planet Nine at ~700 AU has V ≈ 21.7
        // (the prompt's "~22"). Our reflected-light model reuses the shared
        // Neptune-anchored photometry; we pin within ~1 mag of the published
        // value (the V-band reflected model is simplified — Neptune albedo,
        // no atmospheric band structure).
        let v = reflected_v_magnitude(10.0, 700.0);
        assert!(
            (v - reference::V_10ME_700AU).abs() < 1.0,
            "V(10 M⊕, 700 AU) = {v:.2}, published {:.2}",
            reference::V_10ME_700AU
        );
    }

    #[test]
    fn v_brightens_with_mass() {
        // More massive -> larger radius -> more reflected flux -> brighter.
        let v5 = reflected_v_magnitude(5.0, 700.0);
        let v10 = reflected_v_magnitude(10.0, 700.0);
        let v20 = reflected_v_magnitude(20.0, 700.0);
        assert!(
            v20 < v10 && v10 < v5,
            "V: 5={v5:.2} 10={v10:.2} 20={v20:.2}"
        );
    }

    #[test]
    fn v_brightens_at_smaller_distance() {
        let near = reflected_v_magnitude(10.0, 400.0);
        let far = reflected_v_magnitude(10.0, 700.0);
        assert!(
            near < far,
            "near {near:.2} should be brighter than far {far:.2}"
        );
    }

    #[test]
    fn v_dims_by_three_mag_per_doubling() {
        // Reflected light: ~1/(r²Δ²) ≈ 1/r⁴ at large r, so doubling distance
        // dims by ~10 log10(2) ≈ 3.01 mag (the workspace photometry law).
        let v500 = reflected_v_magnitude(10.0, 500.0);
        let v1000 = reflected_v_magnitude(10.0, 1000.0);
        assert_relative_eq!(v1000 - v500, 3.01, epsilon = 0.05);
    }

    #[test]
    fn faint_lowmass_p9_exceeds_optical_survey_limits() {
        // A fainter (lower-mass, more-distant) Planet Nine quickly drops below
        // optical survey limits. PS1 3pi reaches r ≈ 21.5; a 5 M⊕ body at
        // 1000 AU is far fainter in reflected light.
        let ps1 = limiting_magnitude("PS1 3pi").unwrap();
        let v = reflected_v_magnitude(5.0, 1000.0);
        assert!(
            !is_detectable(v, ps1),
            "V(5 M⊕, 1000 AU) = {v:.1} should exceed PS1 limit {ps1}"
        );
        // ... while the nominal 10 M⊕ at 700 AU sits near the deepest optical
        // limit (DES r ≈ 23.8), illustrating how marginal optical detection is.
        let des = limiting_magnitude("DES").unwrap();
        let v_nom = reflected_v_magnitude(10.0, 700.0);
        assert!(
            is_detectable(v_nom, des),
            "V(10 M⊕, 700 AU) = {v_nom:.1} should be within DES limit {des}"
        );
    }

    #[test]
    fn magnitude_table_is_monotone_in_mass() {
        let rows = magnitude_table(&[5.0, 10.0, 20.0], 700.0);
        assert_eq!(rows.len(), 3);
        for w in rows.windows(2) {
            // Heavier -> brighter (smaller mag) in both channels.
            assert!(w[1].v_mag < w[0].v_mag, "V not monotone");
            assert!(w[1].mid_ir_mag < w[0].mid_ir_mag, "mid-IR not monotone");
        }
    }

    #[test]
    fn mid_ir_is_more_favorable_than_optical_for_cold_body() {
        // The paper's central detection message: for the cold (~40-50 K) body
        // the thermal mid-IR (WISE W4, 22 µm) is brighter than reflected
        // optical V, so longward of ~13 µm is the favorable detection regime.
        let row = magnitude_row(10.0, 700.0);
        assert!(
            row.mid_ir_mag < row.v_mag,
            "mid-IR W4 {:.1} should be brighter than V {:.1}",
            row.mid_ir_mag,
            row.v_mag
        );
    }
}
