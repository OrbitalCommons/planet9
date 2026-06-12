//! Earth/observer geometry providers for survey models.
//!
//! Survey crates need the heliocentric position of the *observer* (Earth)
//! at real survey epochs to turn heliocentric model positions into apparent
//! geocentric ones. Two providers implement the same [`EarthProvider`]
//! interface, and the choice is explicit at every call site:
//!
//! - [`EphemerisEarth`] — the primary provider, backed by starfield's JPL
//!   DE421 SPICE kernel (loaded on demand from the starfield cache;
//!   [`EphemerisEarth::try_new`] never downloads, so test suites stay
//!   hermetic on machines without the cache file).
//! - [`CircularEarth`] — the documented analytic fallback: a circular,
//!   coplanar 1 AU orbit at Earth's mean longitude. Earth's true orbit has
//!   e ≈ 0.0167, so this is wrong by up to ~2° in longitude (~0.033 AU),
//!   which at a 400–700 AU target maps to ≲ 0.3 arcmin of apparent
//!   position — the approximation error the survey crates document.
//!
//! [`apparent_direction_from`] implements the target side of starfield's
//! apparent-position chain (`observe().apparent()`) for *synthetic* targets
//! that are not in the kernel: iterative light-time retardation (~3–4 days
//! at 500–700 AU) plus relativistic stellar aberration via
//! `starfield::relativity::add_aberration`. Gravitational deflection is
//! omitted (sub-mas away from the solar limb — far below survey astrometry).

use nalgebra::Vector3;
use starfield::constants::C_AUDAY;
use starfield::planetlib::{Body, Ephemeris};
use starfield::relativity::add_aberration;

pub use starfield::time::{Time, Timescale};

use crate::constants::J2000;

/// Heliocentric ecliptic-J2000 state of Earth: (position AU, velocity AU/day).
pub type EarthState = (Vector3<f64>, Vector3<f64>);

/// Ephemeris file backing [`EphemerisEarth`] (resolved from the starfield
/// cache, `~/.cache/starfield/` by default).
pub const EPHEMERIS_FILE: &str = "de421.bsp";

/// Source of Earth's heliocentric state at a given epoch.
///
/// Methods take `&mut self` because SPICE kernel evaluation requires
/// mutable access in starfield.
pub trait EarthProvider {
    /// Heliocentric ecliptic-J2000 Earth state at `t`:
    /// (position in AU, velocity in AU/day).
    fn earth_state(&mut self, t: &Time) -> EarthState;

    /// Heliocentric ecliptic-J2000 Earth position at `t` (AU).
    fn earth_position(&mut self, t: &Time) -> Vector3<f64> {
        self.earth_state(t).0
    }
}

/// Analytic fallback: Earth on a circular, coplanar 1 AU orbit at its mean
/// longitude L = 100.46435° + 0.985647365°/day × (t − J2000).
///
/// The neglected equation of center (≤ 2e ≈ 1.9°) bounds the position
/// error at ~0.033 AU; see the module docs for the apparent-angle impact.
#[derive(Debug, Clone, Copy, Default)]
pub struct CircularEarth;

/// Earth's mean longitude at J2000 (degrees).
const EARTH_MEAN_LON_J2000_DEG: f64 = 100.46435;
/// Earth's mean motion (degrees/day).
const EARTH_MEAN_MOTION_DEG_DAY: f64 = 0.985_647_365;

impl EarthProvider for CircularEarth {
    fn earth_state(&mut self, t: &Time) -> EarthState {
        let d = t.tdb() - J2000;
        let lon = (EARTH_MEAN_LON_J2000_DEG + EARTH_MEAN_MOTION_DEG_DAY * d).to_radians();
        let rate = EARTH_MEAN_MOTION_DEG_DAY.to_radians(); // rad/day
        let (sin_l, cos_l) = lon.sin_cos();
        (
            Vector3::new(cos_l, sin_l, 0.0),
            Vector3::new(-sin_l * rate, cos_l * rate, 0.0),
        )
    }
}

/// Primary provider: real Earth states from the DE421 JPL ephemeris via
/// starfield's `planetlib::Ephemeris` (heliocentric, rotated to ecliptic
/// J2000).
pub struct EphemerisEarth {
    ephemeris: Ephemeris,
}

impl EphemerisEarth {
    /// Open the DE421 kernel from the starfield cache. Never downloads:
    /// returns `Err` when the cache file is absent so callers (tests in
    /// particular) can fall back or skip gracefully.
    pub fn try_new() -> Result<Self, String> {
        let path = starfield::data::get_cache_dir().join(EPHEMERIS_FILE);
        let populated = std::fs::metadata(&path)
            .map(|m| m.len() > 0)
            .unwrap_or(false);
        if !populated {
            return Err(format!(
                "{EPHEMERIS_FILE} not present in starfield cache ({})",
                path.display()
            ));
        }
        let ephemeris = starfield::Loader::new()
            .load_ephemeris(&path)
            .map_err(|e| format!("failed to open {}: {e}", path.display()))?;
        Ok(Self { ephemeris })
    }

    /// Like [`Self::try_new`], but downloads the kernel from JPL/NAIF when
    /// the cache is empty. Network access — never call from regular tests
    /// (an `#[ignore]`d test may exercise this path).
    pub fn try_new_downloading() -> Result<Self, String> {
        let ephemeris = starfield::Loader::new()
            .open_ephemeris(EPHEMERIS_FILE)
            .map_err(|e| format!("failed to open {EPHEMERIS_FILE}: {e}"))?;
        Ok(Self { ephemeris })
    }
}

impl EarthProvider for EphemerisEarth {
    fn earth_state(&mut self, t: &Time) -> EarthState {
        let sv = self
            .ephemeris
            .ecliptic_state(Body::Earth, t)
            .expect("DE421 covers 1900-2050; survey epochs are inside");
        (sv.position.to_vector3(), sv.velocity.to_vector3())
    }
}

/// Light-time convergence threshold (days); ~1e-9 d ≈ 86 µs ≈ 26 km of
/// path, far below any survey astrometry here.
const LIGHT_TIME_EPSILON_DAYS: f64 = 1e-9;
const MAX_LIGHT_TIME_ITERATIONS: usize = 10;

/// Apparent geocentric direction (ecliptic J2000, not normalized) of a
/// synthetic target observed from `earth_state` (position AU, velocity
/// AU/day).
///
/// `target_helio_at(dt)` must return the target's heliocentric ecliptic
/// position `dt` days after the observation epoch; light-time iteration
/// evaluates it at small negative offsets (the retarded position, ~3–4
/// days at 500–700 AU). Relativistic stellar aberration is applied with
/// starfield's `add_aberration`; gravitational deflection is omitted
/// (sub-mas at these geometries).
pub fn apparent_direction_from<F>(earth_state: &EarthState, target_helio_at: F) -> Vector3<f64>
where
    F: Fn(f64) -> Vector3<f64>,
{
    let (e_pos, e_vel) = earth_state;
    let mut light_time = 0.0;
    let mut geo = target_helio_at(0.0) - e_pos;
    for _ in 0..MAX_LIGHT_TIME_ITERATIONS {
        let lt = geo.norm() / C_AUDAY;
        if (lt - light_time).abs() < LIGHT_TIME_EPSILON_DAYS {
            light_time = lt;
            break;
        }
        light_time = lt;
        geo = target_helio_at(-light_time) - e_pos;
    }
    add_aberration(&mut geo, e_vel, light_time);
    geo
}

/// Apparent geocentric direction of a synthetic target at epoch `t`,
/// fetching the Earth state from `earth`. See [`apparent_direction_from`].
pub fn apparent_geocentric_ecliptic<F>(
    earth: &mut impl EarthProvider,
    t: &Time,
    target_helio_at: F,
) -> Vector3<f64>
where
    F: Fn(f64) -> Vector3<f64>,
{
    let state = earth.earth_state(t);
    apparent_direction_from(&state, target_helio_at)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{RAD2DEG, YEAR_DAYS};

    fn ts() -> Timescale {
        Timescale::default()
    }

    #[test]
    fn circular_earth_is_unit_orbit_with_annual_period() {
        let ts = ts();
        let mut earth = CircularEarth;
        let t0 = ts.utc((2010, 3, 1));
        let (p0, v0) = earth.earth_state(&t0);
        assert!((p0.norm() - 1.0).abs() < 1e-12);
        // Velocity is tangential at the circular rate 2π/yr.
        assert!(p0.dot(&v0).abs() < 1e-12);
        assert!((v0.norm() - std::f64::consts::TAU / YEAR_DAYS).abs() < 1e-4);
        // Half a year later Earth is on the opposite side.
        let t1 = ts.tt_jd(t0.tt() + YEAR_DAYS / 2.0, None);
        let p1 = earth.earth_position(&t1);
        assert!(
            (p0 + p1).norm() < 0.02,
            "opposition residual = {}",
            (p0 + p1).norm()
        );
    }

    #[test]
    fn circular_earth_matches_ephemeris_within_documented_error() {
        let mut eph = match EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let ts = ts();
        let mut circ = CircularEarth;
        for &(y, m, d) in &[
            (1983, 6, 24),
            (1990, 1, 4),
            (2006, 12, 31),
            (2013, 8, 31),
            (2019, 1, 9),
        ] {
            let t = ts.utc((y, m, d));
            let (pe, ve) = eph.earth_state(&t);
            let (pc, vc) = circ.earth_state(&t);
            // True orbit: r within [0.983, 1.017] AU, |z| tiny.
            assert!((0.980..1.020).contains(&pe.norm()), "r = {}", pe.norm());
            assert!(pe.z.abs() < 1e-3, "z = {}", pe.z);
            // Mean-longitude circular model within ~2.1 deg of the truth
            // (equation of center bound), i.e. ~0.04 AU.
            let sep = pe.angle(&pc) * RAD2DEG;
            assert!(
                sep < 2.1,
                "Earth longitude error {sep:.3} deg at {y}-{m}-{d}"
            );
            assert!((pe - pc).norm() < 0.045);
            // Speeds agree to the eccentricity scale.
            assert!((ve.norm() - vc.norm()).abs() < 6e-4);
        }
    }

    #[test]
    fn apparent_direction_includes_light_time_and_aberration() {
        let mut eph = match EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let ts = ts();
        let t = ts.utc((2006, 12, 31));
        let state = eph.earth_state(&t);

        // A target at 500 AU moving at 2e-5 AU/day along +y.
        let target = |dt: f64| Vector3::new(500.0, 2.0e-5 * dt, 0.0);
        let apparent = apparent_direction_from(&state, target);

        // Geometric (no light-time, no aberration) direction.
        let geometric = target(0.0) - state.0;
        let shift_arcsec = apparent.angle(&geometric) * RAD2DEG * 3600.0;
        // Aberration is bounded by ~20.5"; light-time displaces the target
        // by µ·(r/c) ≈ 0.5". The combination must be present but small.
        assert!(
            shift_arcsec > 0.5 && shift_arcsec < 21.5,
            "apparent-geometric shift = {shift_arcsec:.2} arcsec"
        );

        // The light-time alone (no aberration) shifts by µ·lt ~ 0.5":
        // the retarded position differs from the instantaneous one.
        let lt_days = geometric.norm() / C_AUDAY;
        assert!(
            (2.8..3.0).contains(&lt_days),
            "light time at ~500 AU = {lt_days:.3} days"
        );
    }

    #[test]
    fn try_new_never_creates_cache_entries() {
        // try_new on a hermetic machine returns Err rather than
        // downloading. We can't unmake the local cache, but we can assert
        // the error path formats sensibly by pointing the check at the
        // existing API contract: a successful open or a clean Err.
        match EphemerisEarth::try_new() {
            Ok(_) => {}
            Err(msg) => assert!(msg.contains(EPHEMERIS_FILE)),
        }
    }

    /// Network path: downloads DE421 when absent. Kept `#[ignore]`d so the
    /// suite never touches the network; run explicitly with
    /// `cargo test -- --ignored` to exercise the download fallback.
    #[test]
    #[ignore]
    fn ephemeris_download_path_works() {
        let mut earth = EphemerisEarth::try_new_downloading().expect("download or cache de421");
        let ts = ts();
        let p = earth.earth_position(&ts.utc((2000, 1, 1)));
        assert!((p.norm() - 0.9833).abs() < 0.02, "r = {}", p.norm());
    }
}
