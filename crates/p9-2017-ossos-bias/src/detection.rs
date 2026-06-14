//! Survey detection model for large-a TNOs.
//!
//! Shankman et al. (2017) show that a wide-field survey detects a distant
//! eccentric object only when several conditions line up at once, and that
//! the combination imposes "striking biases" on the recovered sample:
//!
//! 1. **Brightness near perihelion.** Reflected sunlight falls as r⁻⁴, so a
//!    150–800 AU object is only above a survey's limiting magnitude for a
//!    small arc near perihelion. We evaluate each object at perihelion
//!    (M = 0, its brightest configuration) using the shared
//!    `p9_core::analysis::photometry` opposition law and a depth from
//!    `p9_core::analysis::surveys`. This is the dominant selection: it ties
//!    detectability to the *perihelion direction*, hence to ϖ.
//!
//! 2. **Ecliptic-latitude / declination acceptance.** Surveys cover a band;
//!    objects whose perihelion sits far off the ecliptic, or below the
//!    survey's declination limit, are missed. We use the ecliptic→equatorial
//!    transform of `p9_core::coords::sky` for the declination cut and a
//!    latitude band for the ecliptic-coverage cut.
//!
//! 3. **Discovery-longitude windows.** A survey only searches certain blocks
//!    of ecliptic longitude (OSSOS's blocks, plus the longitudes other
//!    wide-field surveys happened to cover). An object is discoverable only
//!    if its perihelion ecliptic longitude falls in a covered window. Because
//!    a distant object is bright essentially only at perihelion, the covered
//!    longitudes select a preferred perihelion longitude — and, with the
//!    objects' apsides, a preferred ϖ. This is the mechanism that converts an
//!    isotropic parent into an apparently clustered detected sample.

use p9_core::analysis::photometry::{apparent_magnitude, opposition_delta};
use p9_core::constants::{DEG2RAD, GM_SUN, TWO_PI};
use p9_core::coords::sky::ecliptic_vec_declination;
use p9_core::types::OrbitalElements;

use crate::population::SyntheticEtno;

/// A block of covered ecliptic longitude (radians), with a half-width.
#[derive(Debug, Clone, Copy)]
pub struct LongitudeWindow {
    /// Window centre (ecliptic longitude, radians).
    pub center: f64,
    /// Half-width (radians).
    pub half_width: f64,
}

impl LongitudeWindow {
    fn contains(&self, lon: f64) -> bool {
        let mut d = (lon - self.center).rem_euclid(TWO_PI);
        if d > std::f64::consts::PI {
            d -= TWO_PI;
        }
        d.abs() <= self.half_width
    }
}

/// Parameters of the wide-field survey detection model.
#[derive(Debug, Clone)]
pub struct SurveyModel {
    /// Limiting magnitude (V/r). Default from the deep wide-field surveys
    /// (`p9_core::analysis::surveys` DES depth, the deepest tabulated).
    pub limiting_mag: f64,
    /// Declination floor (radians): objects below this are unobservable from
    /// the (predominantly northern) discovery sites.
    pub dec_min: f64,
    /// Ecliptic-latitude half-band (radians): coverage falls off beyond this.
    pub ecliptic_lat_band: f64,
    /// Covered discovery-longitude windows in ecliptic longitude.
    pub windows: Vec<LongitudeWindow>,
}

impl Default for SurveyModel {
    fn default() -> Self {
        // Depth: use the deepest tabulated wide-field survey (DES, r ≈ 23.8)
        // from the shared survey table.
        let limiting_mag = p9_core::analysis::surveys::limiting_magnitude("DES").unwrap_or(23.8);

        // Discovery-longitude windows. OSSOS observed a handful of blocks at
        // specific ecliptic longitudes; combined with the other wide-field
        // surveys that contributed large-a TNOs, the effective coverage is a
        // few longitude blocks, NOT the full circle, and those blocks are
        // *clumped* on the sky rather than evenly spaced (the deep large-a
        // discoveries came from a limited set of nearby fields). That clumping
        // — not the mere fact of incomplete coverage — is what manufactures a
        // single dominant apparent-ϖ lobe out of an isotropic parent: evenly
        // spaced blocks would seed several antipodal lobes and wash R̄ back
        // toward zero. Centres mimic a clumped block layout around
        // λ ≈ 0°–60° with one offset block, each ~20° half-width.
        let windows = [0.0, 30.0, 60.0, 210.0]
            .iter()
            .map(|&c| LongitudeWindow {
                center: c * DEG2RAD,
                half_width: 20.0 * DEG2RAD,
            })
            .collect();

        Self {
            limiting_mag,
            dec_min: -20.0 * DEG2RAD,
            ecliptic_lat_band: 20.0 * DEG2RAD,
            windows,
        }
    }
}

/// State of an orbit placed exactly at perihelion (M = 0): its brightest,
/// most-detectable configuration.
fn at_perihelion(o: &SyntheticEtno) -> OrbitalElements {
    OrbitalElements {
        mean_anomaly: 0.0,
        ..o.elements
    }
}

/// Why an object was or was not detected (for diagnostics).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectionOutcome {
    Detected,
    TooFaint,
    BelowDeclination,
    OffEclipticBand,
    OutsideLongitudeWindows,
}

/// Decide whether a synthetic object would be detected, evaluating it at
/// perihelion (its brightest point — a survey that covers that longitude
/// would catch it there).
pub fn detection_outcome(o: &SyntheticEtno, model: &SurveyModel) -> DetectionOutcome {
    let elem = at_perihelion(o);
    let state = elem.to_state_vector(GM_SUN);
    let r = state.pos.norm();
    let lon = state.pos.y.atan2(state.pos.x).rem_euclid(TWO_PI);
    let beta = (state.pos.z / r).asin();
    let dec = ecliptic_vec_declination(&state.pos);

    // Brightness: apparent magnitude at opposition (Δ = r − 1 AU).
    let m = apparent_magnitude(o.h_mag, r, opposition_delta(r));
    if m > model.limiting_mag {
        return DetectionOutcome::TooFaint;
    }

    if dec < model.dec_min {
        return DetectionOutcome::BelowDeclination;
    }

    if beta.abs() > model.ecliptic_lat_band {
        return DetectionOutcome::OffEclipticBand;
    }

    if !model.windows.iter().any(|w| w.contains(lon)) {
        return DetectionOutcome::OutsideLongitudeWindows;
    }

    DetectionOutcome::Detected
}

/// True iff the object is detected.
pub fn is_detected(o: &SyntheticEtno, model: &SurveyModel) -> bool {
    detection_outcome(o, model) == DetectionOutcome::Detected
}

/// Filter a population to the detected subsample.
pub fn detected_subsample(pop: &[SyntheticEtno], model: &SurveyModel) -> Vec<SyntheticEtno> {
    pop.iter()
        .copied()
        .filter(|o| is_detected(o, model))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::population::{generate_population, PopulationParams};
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_faint_distant_object_rejected() {
        // Faint (H = 9) object with a far perihelion: too faint at perihelion.
        let o = SyntheticEtno {
            elements: OrbitalElements {
                a: 800.0,
                e: 0.9,
                i: 0.1,
                omega_big: 0.2,
                omega: 0.1,
                mean_anomaly: 0.0,
            },
            h_mag: 9.0,
        };
        // q = 80 AU, H = 9 → m well past DES depth.
        assert_eq!(
            detection_outcome(&o, &SurveyModel::default()),
            DetectionOutcome::TooFaint
        );
    }

    #[test]
    fn test_some_objects_are_detected() {
        let mut rng = StdRng::seed_from_u64(3);
        let pop = generate_population(20_000, &PopulationParams::default(), &mut rng);
        let det = detected_subsample(&pop, &SurveyModel::default());
        assert!(
            det.len() > 20 && det.len() < pop.len(),
            "detected {} of {} — selection should be strong but non-empty",
            det.len(),
            pop.len()
        );
    }

    #[test]
    fn test_longitude_window_contains_wraps() {
        let w = LongitudeWindow {
            center: 5.0 * DEG2RAD,
            half_width: 10.0 * DEG2RAD,
        };
        assert!(w.contains(358.0 * DEG2RAD));
        assert!(w.contains(12.0 * DEG2RAD));
        assert!(!w.contains(30.0 * DEG2RAD));
    }
}
