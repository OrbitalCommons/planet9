//! Giant planet initial conditions.
//!
//! Two sources, one [`MassiveBody`] convention (heliocentric ecliptic J2000,
//! AU, AU/day):
//!
//! - Hard-coded J2000.0 states ([`giant_planets_j2000`] and friends) — the
//!   hermetic default for tests and reproducible runs. The values are JPL
//!   DE440 planet-system *barycenter* states evaluated at J2000.0
//!   (JD 2451545.0 TDB), pinned against the live ephemeris by a kernel-gated
//!   test below.
//! - Ephemeris states at arbitrary epochs ([`PlanetEphemeris`],
//!   [`giant_planets_at`]) — backed by starfield's SPK reader over the DE
//!   kernels in the starfield cache (`~/.cache/starfield/`). Loading never
//!   downloads (see [`PlanetEphemeris::try_new`]); an explicit
//!   [`PlanetEphemeris::try_new_downloading`] exists for non-hermetic use.
//!
//! Barycenter vs body center: starfield's `planetlib` resolves the outer
//! planets to their NAIF *system barycenters* ("jupiter barycenter", …,
//! "neptune barycenter") — the standard targets in the planetary DE kernels.
//! That is the right object for our dynamics: p9-core's `GM_*` /
//! `MASS_*_SOLAR` constants are DE440 planetary-*system* GMs (planet +
//! satellites), so barycenter states and system masses describe the same
//! body. The center-vs-barycenter offset is tiny anyway (≲ 1e-6 AU,
//! ≲ 3e-7 AU/day even for Jupiter/Neptune).
//!
//! Masses vs states: the ephemeris provides *states only*. GM and mass
//! always come from p9-core's constants (used by the integrators), never
//! from the kernel's own GM header — so hermetic (J2000) and
//! ephemeris-seeded integrations share one mass model and differ only in
//! the starting state.

use nalgebra::Vector3;
use starfield::planetlib::{Body, Ephemeris};

pub use starfield::time::{Time, Timescale};

use crate::constants::*;
use crate::types::{MassiveBody, StateVector};

/// DE kernels usable for planet states, in preference order, resolved from
/// the starfield cache (`~/.cache/starfield/` by default). DE440 covers
/// 1550–2650; DE421 (also used by `coords::observer`) covers 1900–2050.
pub const PLANET_EPHEMERIS_FILES: [&str; 2] = ["de440.bsp", "de421.bsp"];

/// Attach p9-core's physical parameters (system GM/mass, radius) to a
/// giant-planet state. Single source for both the hard-coded J2000 states
/// and the ephemeris path, so the two can never disagree on masses.
/// (Planetary J2/J4 enter the dynamics only through the solar-quadrupole
/// proxy in `crate::forces::j2_secular`, not per-body fields.)
fn giant_body(body: Body, state: StateVector) -> MassiveBody {
    let (name, gm, mass, radius_au) = match body {
        Body::Jupiter => ("Jupiter", GM_JUPITER, MASS_JUPITER_SOLAR, RADIUS_JUPITER_AU),
        Body::Saturn => ("Saturn", GM_SATURN, MASS_SATURN_SOLAR, RADIUS_SATURN_AU),
        Body::Uranus => ("Uranus", GM_URANUS, MASS_URANUS_SOLAR, RADIUS_URANUS_AU),
        Body::Neptune => ("Neptune", GM_NEPTUNE, MASS_NEPTUNE_SOLAR, RADIUS_NEPTUNE_AU),
        other => panic!("{} is not a giant planet", other.name()),
    };
    MassiveBody {
        name: name.to_string(),
        gm,
        mass,
        state,
        radius_au,
    }
}

// The J2000 constants below are JPL DE440 planet-system barycenter states
// relative to the Sun at JD 2451545.0 TDB, rotated to ecliptic J2000
// (extracted via starfield's SPK reader from de440.bsp and pinned by the
// kernel-gated `j2000_constants_match_ephemeris` test). The previous values
// (pre issue #40) disagreed with DE440 by ~7e-3 AU / ~1e-5 AU/day — they
// had evidently been transcribed from an approximate source.

/// Jupiter state vector at J2000.0, heliocentric ecliptic J2000 (DE440).
pub fn jupiter_j2000() -> MassiveBody {
    giant_body(
        Body::Jupiter,
        StateVector::new(
            Vector3::new(
                4.001_177_161_126_056_7,
                2.938_576_106_966_103,
                -0.101_785_276_633_744_04,
            ),
            Vector3::new(
                -0.004_568_313_526_752_718,
                0.006_443_206_010_221_306,
                0.000_075_579_520_257_744_97,
            ),
        ),
    )
}

/// Saturn state vector at J2000.0, heliocentric ecliptic J2000 (DE440).
pub fn saturn_j2000() -> MassiveBody {
    giant_body(
        Body::Saturn,
        StateVector::new(
            Vector3::new(
                6.406_408_859_532_808,
                6.569_989_616_912_556,
                -0.369_076_426_244_853_4,
            ),
            Vector3::new(
                -0.004_292_351_873_843_249,
                0.003_890_315_698_014_864,
                0.000_102_947_849_398_980_56,
            ),
        ),
    )
}

/// Uranus state vector at J2000.0, heliocentric ecliptic J2000 (DE440).
pub fn uranus_j2000() -> MassiveBody {
    giant_body(
        Body::Uranus,
        StateVector::new(
            Vector3::new(
                14.431_856_614_381_783,
                -13.734_321_468_770_963,
                -0.238_141_687_799_250_16,
            ),
            Vector3::new(
                0.002_678_105_591_420_150_3,
                0.002_672_695_895_712_686_4,
                -0.000_024_772_063_555_051_14,
            ),
        ),
    )
}

/// Neptune state vector at J2000.0, heliocentric ecliptic J2000 (DE440).
pub fn neptune_j2000() -> MassiveBody {
    giant_body(
        Body::Neptune,
        StateVector::new(
            Vector3::new(
                16.812_046_968_052_883,
                -24.991_762_889_284_69,
                0.127_222_879_920_248_39,
            ),
            Vector3::new(
                0.002_579_274_259_047_737,
                0.001_776_900_820_016_125_8,
                -0.000_095_909_614_985_915_96,
            ),
        ),
    )
}

/// All four giant planets at J2000.0.
pub fn giant_planets_j2000() -> Vec<MassiveBody> {
    vec![
        jupiter_j2000(),
        saturn_j2000(),
        uranus_j2000(),
        neptune_j2000(),
    ]
}

/// Neptune only (many papers use J2 for inner giants + direct Neptune).
pub fn neptune_only_j2000() -> Vec<MassiveBody> {
    vec![neptune_j2000()]
}

/// Giant-planet states at arbitrary epochs from a JPL DE kernel.
///
/// Methods take `&mut self` because SPICE kernel evaluation requires
/// mutable access in starfield. See the module docs for the
/// barycenter/mass conventions.
pub struct PlanetEphemeris {
    ephemeris: Ephemeris,
}

impl PlanetEphemeris {
    /// Open the first available DE kernel from the starfield cache
    /// ([`PLANET_EPHEMERIS_FILES`] preference order). Never downloads:
    /// returns `Err` when no kernel is cached so callers (tests in
    /// particular) can fall back to [`giant_planets_j2000`] or skip.
    pub fn try_new() -> Result<Self, String> {
        let cache = starfield::data::get_cache_dir();
        for file in PLANET_EPHEMERIS_FILES {
            let path = cache.join(file);
            let populated = std::fs::metadata(&path)
                .map(|m| m.len() > 0)
                .unwrap_or(false);
            if !populated {
                continue;
            }
            let ephemeris = starfield::Loader::new()
                .load_ephemeris(&path)
                .map_err(|e| format!("failed to open {}: {e}", path.display()))?;
            return Ok(Self { ephemeris });
        }
        Err(format!(
            "no DE kernel ({}) in starfield cache ({})",
            PLANET_EPHEMERIS_FILES.join(", "),
            cache.display()
        ))
    }

    /// Like [`Self::try_new`], but downloads DE440 from JPL/NAIF when the
    /// cache is empty. Network access — never call from regular tests (an
    /// `#[ignore]`d test may exercise this path).
    pub fn try_new_downloading() -> Result<Self, String> {
        let file = PLANET_EPHEMERIS_FILES[0];
        let ephemeris = starfield::Loader::new()
            .open_ephemeris(file)
            .map_err(|e| format!("failed to open {file}: {e}"))?;
        Ok(Self { ephemeris })
    }

    /// One giant planet at `t`: heliocentric ecliptic-J2000 barycenter state
    /// from the kernel, physical parameters from p9-core constants.
    pub fn giant_planet_at(&mut self, body: Body, t: &Time) -> Result<MassiveBody, String> {
        let sv = self
            .ephemeris
            .ecliptic_state(body, t)
            .map_err(|e| format!("ephemeris lookup for {} failed: {e}", body.name()))?;
        Ok(giant_body(
            body,
            StateVector::new(sv.position.to_vector3(), sv.velocity.to_vector3()),
        ))
    }

    /// All four giant planets at `t` (ephemeris analogue of
    /// [`giant_planets_j2000`]).
    pub fn giant_planets_at(&mut self, t: &Time) -> Result<Vec<MassiveBody>, String> {
        [Body::Jupiter, Body::Saturn, Body::Uranus, Body::Neptune]
            .into_iter()
            .map(|b| self.giant_planet_at(b, t))
            .collect()
    }

    /// Neptune only at `t` (ephemeris analogue of [`neptune_only_j2000`]).
    pub fn neptune_only_at(&mut self, t: &Time) -> Result<Vec<MassiveBody>, String> {
        Ok(vec![self.giant_planet_at(Body::Neptune, t)?])
    }
}

/// All four giant planets at epoch `t` from the cached DE kernel.
///
/// One-shot convenience over [`PlanetEphemeris`] (opens the kernel per
/// call — hold a `PlanetEphemeris` when querying many epochs). Never
/// downloads; `Err` when no kernel is cached.
pub fn giant_planets_at(t: &Time) -> Result<Vec<MassiveBody>, String> {
    PlanetEphemeris::try_new()?.giant_planets_at(t)
}

/// Neptune only at epoch `t` from the cached DE kernel. See
/// [`giant_planets_at`].
pub fn neptune_only_at(t: &Time) -> Result<Vec<MassiveBody>, String> {
    PlanetEphemeris::try_new()?.neptune_only_at(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrator::whm::{barycentric_energy, WhmIntegrator};
    use crate::types::SimConfig;
    use uom::si::length::astronomical_unit;

    fn ephemeris_or_skip() -> Option<PlanetEphemeris> {
        match PlanetEphemeris::try_new() {
            Ok(e) => Some(e),
            Err(_) => {
                eprintln!("skipping: no DE kernel in starfield cache");
                None
            }
        }
    }

    #[test]
    fn hermetic_j2000_states_are_sane() {
        // No kernel needed: distances and speeds near the catalog values.
        let expected_au = [
            (5.2, "Jupiter"),
            (9.5, "Saturn"),
            (19.2, "Uranus"),
            (30.1, "Neptune"),
        ];
        for (body, (a, name)) in giant_planets_j2000().iter().zip(expected_au) {
            assert_eq!(body.name, name);
            let r = body.state.distance_typed().get::<astronomical_unit>();
            assert!(
                (r - a).abs() / a < 0.06,
                "{name} at r = {r:.2} AU, expected near {a} AU"
            );
            let vis_viva = (GM_SUN * (2.0 / r - 1.0 / a)).sqrt();
            assert!((body.state.speed() - vis_viva).abs() / vis_viva < 0.05);
        }
    }

    /// Pin the hard-coded J2000 constants against the live DE kernel at
    /// JD 2451545.0 TDB (kernel-gated).
    #[test]
    fn j2000_constants_match_ephemeris() {
        let Some(mut eph) = ephemeris_or_skip() else {
            return;
        };
        let ts = Timescale::default();
        let t = ts.tdb_jd(J2000);
        let live = eph.giant_planets_at(&t).expect("J2000 inside DE coverage");

        for (pinned, live) in giant_planets_j2000().iter().zip(&live) {
            let dp = (pinned.state.pos - live.state.pos).norm();
            let dv = (pinned.state.vel - live.state.vel).norm();
            eprintln!(
                "{}: |dpos| = {dp:.3e} AU, |dvel| = {dv:.3e} AU/day",
                pinned.name
            );
            // The constants are DE440 barycenter states; allow for DE421
            // (the fallback kernel) differing from DE440 by up to a few
            // hundred km (~5e-6 AU) for the outer planets.
            assert!(
                dp < 5e-6,
                "{}: |dpos| = {dp:.3e} AU vs ephemeris at J2000",
                pinned.name
            );
            assert!(
                dv < 5e-9,
                "{}: |dvel| = {dv:.3e} AU/day vs ephemeris at J2000",
                pinned.name
            );
        }
    }

    /// Epoch-flexible start: integrate the four giants from a 2017-01-01
    /// ephemeris state for 10 kyr and require energy conservation
    /// comparable to the hermetic J2000 start (kernel-gated).
    #[test]
    fn epoch_start_conserves_energy_like_j2000_start() {
        if ephemeris_or_skip().is_none() {
            return;
        }
        let ts = Timescale::default();
        let epoch = ts.utc((2017, 1, 1));

        let max_de = |mut bodies: Vec<MassiveBody>| -> f64 {
            let config = SimConfig {
                dt: 200.0, // < P_Jupiter / 20
                t_start: 0.0,
                t_end: 10_000.0 * YEAR_DAYS,
                removal_inner_au: 0.1,
                removal_outer_au: 1e6,
                snapshot_interval_days: f64::INFINITY,
                hybrid_changeover_hill: 3.0,
                bs_epsilon: 1e-12,
            };
            let integrator = WhmIntegrator::new();
            let mut particles = vec![];
            let mut active = vec![];
            let e0 = barycentric_energy(&bodies);
            let n_steps = (config.t_end / config.dt).round() as usize;
            let mut worst: f64 = 0.0;
            for step in 0..n_steps {
                integrator.step(&mut bodies, &mut particles, &mut active, config.dt, &config);
                if step % 500 == 499 {
                    let e1 = barycentric_energy(&bodies);
                    worst = worst.max(((e1 - e0) / e0).abs());
                }
            }
            worst
        };

        let (_, epoch_bodies) = SimConfig::standard_4gyr()
            .with_bodies_at_epoch(&epoch)
            .expect("kernel present (checked above)");
        let de_epoch = max_de(epoch_bodies);
        let de_j2000 = max_de(giant_planets_j2000());

        assert!(de_j2000 < 1e-5, "J2000 start dE/E = {de_j2000:.2e}");
        assert!(de_epoch < 1e-5, "2017 start dE/E = {de_epoch:.2e}");
        // "Comparable": same symplectic error band, within an order of
        // magnitude of each other.
        let ratio = de_epoch / de_j2000;
        assert!(
            (0.1..10.0).contains(&ratio),
            "energy-error ratio epoch/J2000 = {ratio:.2}"
        );
    }

    #[test]
    fn neptune_only_at_matches_giant_planets_at() {
        let Some(mut eph) = ephemeris_or_skip() else {
            return;
        };
        let ts = Timescale::default();
        let t = ts.utc((2010, 6, 1));
        let nep = &eph.neptune_only_at(&t).unwrap()[0];
        let all = eph.giant_planets_at(&t).unwrap();
        assert_eq!(nep.name, "Neptune");
        assert_eq!(all.len(), 4);
        assert!((nep.state.pos - all[3].state.pos).norm() < 1e-12);
        // Physical parameters come from p9-core constants, not the kernel.
        assert_eq!(nep.gm, GM_NEPTUNE);
        assert_eq!(nep.mass, MASS_NEPTUNE_SOLAR);
    }

    #[test]
    fn try_new_never_creates_cache_entries() {
        // Hermetic contract: a successful open or a clean Err naming the
        // kernels — never a download.
        match PlanetEphemeris::try_new() {
            Ok(_) => {}
            Err(msg) => assert!(msg.contains(PLANET_EPHEMERIS_FILES[0])),
        }
    }

    /// Network path: downloads DE440 when absent. Kept `#[ignore]`d so the
    /// suite never touches the network; run explicitly with
    /// `cargo test -- --ignored`.
    #[test]
    #[ignore]
    fn ephemeris_download_path_works() {
        let mut eph = PlanetEphemeris::try_new_downloading().expect("download or cache de440");
        let ts = Timescale::default();
        let bodies = eph.giant_planets_at(&ts.tdb_jd(J2000)).unwrap();
        assert!((bodies[0].state.distance_typed().get::<astronomical_unit>() - 5.0).abs() < 0.5);
    }
}
