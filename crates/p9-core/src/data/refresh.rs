//! Live JPL SBDB refresh / drift detection for the frozen element tables.
//!
//! The workspace's observational tables (`data::etno::BROWN_2017_SAMPLE`, the
//! p9-2024-neptune-crossing sample) are *frozen snapshots*: provenance-commented
//! constants that regression tests pin against, so default builds and tests are
//! deterministic and offline. They were originally transcribed by hand, and
//! hand transcription is exactly how a ~30 m near-Earth asteroid once ended up
//! listed as a 110 AU TNO (MODELING_REVIEW.md §5 N1).
//!
//! This module closes that hole: it diffs a snapshot table against live SBDB
//! lookups, element by element, with per-element tolerances. The pure diff
//! logic is always compiled (and unit-tested offline); the functions that
//! actually talk to JPL are gated behind the `sbdb-refresh` cargo feature and
//! exercised only by `#[ignore]`d integration tests.
//!
//! # Tolerances
//!
//! Orbit solutions drift between fits: as the observed arc lengthens JPL
//! republishes elements at new epochs, and for the distant high-eccentricity
//! objects in these tables the perihelion distance and the angles are well
//! constrained while `a` and `e` are strongly correlated and poorly
//! constrained — fractional shifts in `a` of a few ×0.1% (occasionally
//! percents for the longest-period objects, e.g. Sedna) and fractions of a
//! degree in the angles are normal between solution epochs. A tolerance
//! therefore has two parts: an allowance for solution drift plus the rounding
//! half-width of the snapshot itself. See [`Tolerances::brown2017`] and
//! [`Tolerances::full_precision`] for the two concrete choices and their
//! rationale.
//!
//! A diff that exceeds tolerance means one of two things, both of which a
//! human should look at: the snapshot was mistranscribed (fix it), or the JPL
//! solution genuinely moved (update the snapshot and its provenance comment,
//! re-pinning any downstream regression values that depend on it).

use crate::constants::DEG2RAD;
use crate::data::etno::BROWN_2017_SAMPLE;
use crate::types::solve_kepler;

/// Which orbital element a diff refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Element {
    /// Semi-major axis a (AU)
    SemiMajorAxis,
    /// Eccentricity e
    Eccentricity,
    /// Inclination i (degrees)
    Inclination,
    /// Argument of perihelion ω (degrees)
    ArgPerihelion,
    /// Longitude of ascending node Ω (degrees)
    AscendingNode,
    /// Absolute magnitude H
    AbsoluteMagnitude,
}

impl Element {
    /// Short human-readable label.
    pub fn label(&self) -> &'static str {
        match self {
            Element::SemiMajorAxis => "a [AU]",
            Element::Eccentricity => "e",
            Element::Inclination => "i [deg]",
            Element::ArgPerihelion => "omega [deg]",
            Element::AscendingNode => "Omega [deg]",
            Element::AbsoluteMagnitude => "H [mag]",
        }
    }

    /// Whether the element is an angle in degrees (delta is wrap-aware).
    pub fn is_angle(&self) -> bool {
        matches!(
            self,
            Element::Inclination | Element::ArgPerihelion | Element::AscendingNode
        )
    }
}

/// One row of a frozen table, in the form the differ compares. Elements the
/// table does not carry are `None` and are skipped.
#[derive(Debug, Clone, Copy)]
pub struct ElementSnapshot {
    /// Human-readable name used in diff reports.
    pub name: &'static str,
    /// SBDB search string (`sstr`) used for the live lookup.
    pub designation: &'static str,
    /// Semi-major axis (AU)
    pub a: Option<f64>,
    /// Eccentricity
    pub e: Option<f64>,
    /// Inclination (degrees)
    pub i_deg: Option<f64>,
    /// Argument of perihelion (degrees)
    pub omega_deg: Option<f64>,
    /// Longitude of ascending node (degrees)
    pub omega_big_deg: Option<f64>,
    /// Absolute magnitude H
    pub h_mag: Option<f64>,
}

/// Osculating elements (and the orbit-propagation extras) from a live SBDB
/// lookup, decoupled from the starfield response types so the diff logic can
/// be unit-tested offline.
#[derive(Debug, Clone, Default)]
pub struct LiveElements {
    /// Semi-major axis (AU)
    pub a: Option<f64>,
    /// Eccentricity
    pub e: Option<f64>,
    /// Inclination (degrees)
    pub i_deg: Option<f64>,
    /// Argument of perihelion (degrees)
    pub omega_deg: Option<f64>,
    /// Longitude of ascending node (degrees)
    pub omega_big_deg: Option<f64>,
    /// Absolute magnitude H
    pub h_mag: Option<f64>,
    /// Epoch of osculation (JD TDB)
    pub epoch_jd: Option<f64>,
    /// Mean anomaly at epoch (degrees)
    pub mean_anomaly_deg: Option<f64>,
    /// Mean motion (degrees/day)
    pub mean_motion_deg_per_day: Option<f64>,
    /// Date of the first observation in the fitted arc (SBDB `first_obs`,
    /// "YYYY-MM-DD"); equals the discovery date except where precovery
    /// observations were later linked.
    pub first_obs: Option<String>,
    /// JPL orbit solution ID (provenance for reports)
    pub orbit_id: Option<String>,
}

/// Per-element comparison tolerances. Each is the *total* allowed
/// |live − snapshot|: solution drift allowance plus the snapshot's rounding
/// half-width.
#[derive(Debug, Clone, Copy)]
pub struct Tolerances {
    /// Fractional tolerance on a (relative to the snapshot value).
    pub a_frac: f64,
    /// Absolute floor on the a tolerance (AU) — covers table rounding.
    pub a_abs: f64,
    /// Absolute tolerance on e.
    pub e_abs: f64,
    /// Absolute tolerance on i, ω, Ω (degrees, wrap-aware).
    pub angle_deg: f64,
    /// Absolute tolerance on H (mag).
    pub h_mag: f64,
}

impl Tolerances {
    /// Tolerances for `data::etno::BROWN_2017_SAMPLE`.
    ///
    /// That table rounds a to the AU (±0.5), e to 0.01 (±0.005), angles to
    /// 0.1° (±0.05), H to 0.1 (±0.05). On top of the rounding:
    /// - a: 1% fractional drift. For these a ≳ 230 AU, e ≳ 0.8 orbits a is
    ///   the loosest-constrained element (strongly correlated with e); a few
    ///   ×0.1% is routine between solution epochs and ~1% is reached by the
    ///   longest-period members. Anything beyond 1% is worth a human look.
    /// - e: 0.01 — drift in e tracks a/q ≈ (1−e) at the few ×0.001 level,
    ///   comparable to the table's own rounding.
    /// - angles: 0.5° — i, ω, Ω are typically stable to ≲0.1° once an orbit
    ///   is multi-opposition; 0.5° flags real solution changes, not noise.
    /// - H: 0.5 mag — SBDB absolute magnitudes are re-derived photometric
    ///   fits and routinely move by a few tenths, independent of astrometry.
    pub fn brown2017() -> Self {
        Tolerances {
            a_frac: 0.01,
            a_abs: 1.0,
            e_abs: 0.01,
            angle_deg: 0.5,
            h_mag: 0.5,
        }
    }

    /// Tolerances for tables transcribed at (near) full SBDB precision, such
    /// as the p9-2024-neptune-crossing sample (a to 0.1 AU, e to 1e-4, i to
    /// 0.01°). Drift allowances are the same physics as
    /// [`Tolerances::brown2017`] but tighter because the snapshot is recent
    /// (June 2026) and unrounded: 0.5% in a (these orbits reach a ~ 800 AU
    /// with short arcs), 0.005 in e, 0.1° in angles.
    pub fn full_precision() -> Self {
        Tolerances {
            a_frac: 0.005,
            a_abs: 0.2,
            e_abs: 0.005,
            angle_deg: 0.1,
            h_mag: 0.5,
        }
    }

    fn for_element(&self, element: Element, snapshot_value: f64) -> f64 {
        match element {
            Element::SemiMajorAxis => (self.a_frac * snapshot_value.abs()).max(self.a_abs),
            Element::Eccentricity => self.e_abs,
            Element::Inclination | Element::ArgPerihelion | Element::AscendingNode => {
                self.angle_deg
            }
            Element::AbsoluteMagnitude => self.h_mag,
        }
    }
}

/// One out-of-tolerance element: which object, which element, what the frozen
/// table says, what JPL says now, and by how much they disagree.
#[derive(Debug, Clone)]
pub struct EtnoDiff {
    /// Object name (from the snapshot table).
    pub object: String,
    /// Which element drifted.
    pub element: Element,
    /// Value in the frozen table.
    pub snapshot: f64,
    /// Live SBDB value.
    pub live: f64,
    /// live − snapshot (wrap-aware for angles, in the element's own units).
    pub delta: f64,
    /// The tolerance that was exceeded (same units).
    pub tolerance: f64,
    /// JPL orbit solution ID of the live lookup, if reported.
    pub orbit_id: Option<String>,
}

impl std::fmt::Display for EtnoDiff {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: {} snapshot {} vs live {} (delta {:+.4}, tolerance {:.4}{})",
            self.object,
            self.element.label(),
            self.snapshot,
            self.live,
            self.delta,
            self.tolerance,
            match &self.orbit_id {
                Some(id) => format!(", JPL soln {id}"),
                None => String::new(),
            }
        )
    }
}

/// Signed angle difference live − snapshot wrapped to (−180°, 180°].
fn angle_delta_deg(snapshot: f64, live: f64) -> f64 {
    let mut d = (live - snapshot).rem_euclid(360.0);
    if d > 180.0 {
        d -= 360.0;
    }
    d
}

/// Diff one snapshot row against live elements. Returns only the elements
/// whose |delta| exceeds tolerance; elements missing on either side are
/// skipped (the snapshot tables only carry `Some` for vetted columns, and
/// SBDB always reports the Keplerian set for real objects).
pub fn diff_elements(
    snapshot: &ElementSnapshot,
    live: &LiveElements,
    tol: &Tolerances,
) -> Vec<EtnoDiff> {
    let pairs = [
        (Element::SemiMajorAxis, snapshot.a, live.a),
        (Element::Eccentricity, snapshot.e, live.e),
        (Element::Inclination, snapshot.i_deg, live.i_deg),
        (Element::ArgPerihelion, snapshot.omega_deg, live.omega_deg),
        (
            Element::AscendingNode,
            snapshot.omega_big_deg,
            live.omega_big_deg,
        ),
        (Element::AbsoluteMagnitude, snapshot.h_mag, live.h_mag),
    ];

    let mut diffs = Vec::new();
    for (element, snap_val, live_val) in pairs {
        let (Some(s), Some(l)) = (snap_val, live_val) else {
            continue;
        };
        let delta = if element.is_angle() {
            angle_delta_deg(s, l)
        } else {
            l - s
        };
        let tolerance = tol.for_element(element, s);
        if delta.abs() > tolerance {
            diffs.push(EtnoDiff {
                object: snapshot.name.to_string(),
                element,
                snapshot: s,
                live: l,
                delta,
                tolerance,
                orbit_id: live.orbit_id.clone(),
            });
        }
    }
    diffs
}

/// The Brown (2017) ETNO table as snapshot rows for the differ.
pub fn etno_snapshots() -> Vec<ElementSnapshot> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| ElementSnapshot {
            name: k.name,
            designation: k.name,
            a: Some(k.a),
            e: Some(k.e),
            i_deg: Some(k.i_deg),
            omega_deg: Some(k.omega_deg),
            omega_big_deg: Some(k.omega_big_deg),
            h_mag: Some(k.h_mag),
        })
        .collect()
}

/// Julian Date (UT) of a calendar date with fractional day (proleptic
/// Gregorian; standard Meeus algorithm). `calendar_to_jd(2000, 1, 1.5)` is
/// J2000.0 = 2451545.0.
pub fn calendar_to_jd(year: i32, month: u32, day: f64) -> f64 {
    let (y, m) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };
    let a = (y as f64 / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (m as f64 + 1.0)).floor() + day + b - 1524.5
}

/// Parse an SBDB date string ("YYYY-MM-DD", optionally with a fractional
/// day) to a Julian Date at 0h UT (plus the fraction if present).
pub fn parse_sbdb_date(s: &str) -> Option<f64> {
    let mut parts = s.trim().splitn(3, '-');
    let year: i32 = parts.next()?.parse().ok()?;
    let month: u32 = parts.next()?.parse().ok()?;
    let day: f64 = parts.next()?.parse().ok()?;
    if !(1..=12).contains(&month) || !(1.0..32.0).contains(&day) {
        return None;
    }
    Some(calendar_to_jd(year, month, day))
}

/// Heliocentric distance (AU) of an object at epoch `jd`, propagated on its
/// own (two-body) live elements: M(t) = M₀ + n·(t − t₀), Kepler-solved to the
/// eccentric anomaly, r = a(1 − e·cos E). Good to the two-body level over the
/// decade-scale spans between a solution epoch and a discovery date, which is
/// far better than the ~15% heuristic it replaces.
pub fn heliocentric_distance_at(live: &LiveElements, jd: f64) -> Option<f64> {
    let a = live.a?;
    let e = live.e?;
    let m0 = live.mean_anomaly_deg?;
    let n = live.mean_motion_deg_per_day?;
    let epoch = live.epoch_jd?;
    if !(0.0..1.0).contains(&e) || a <= 0.0 {
        return None;
    }
    let m = (m0 + n * (jd - epoch)) * DEG2RAD;
    let ea = solve_kepler(e, m);
    Some(a * (1.0 - e * ea.cos()))
}

/// Heliocentric distance (AU) at the first observation of the fitted arc
/// (SBDB `first_obs`) — the live-data replacement for "discovery distance".
/// `None` if the lookup lacked any of the needed fields.
pub fn first_obs_distance(live: &LiveElements) -> Option<f64> {
    let jd = parse_sbdb_date(live.first_obs.as_deref()?)?;
    heliocentric_distance_at(live, jd)
}

/// Live-fetch helpers (network). Gated behind the `sbdb-refresh` feature so
/// default builds and tests never touch the network.
#[cfg(feature = "sbdb-refresh")]
mod live {
    use super::*;
    use starfield::sbdb::SbdbClient;

    /// Look up one object on the live SBDB and flatten the response into
    /// [`LiveElements`].
    pub fn fetch_live_elements(
        client: &SbdbClient,
        designation: &str,
    ) -> Result<LiveElements, String> {
        let resp = client
            .lookup(designation)
            .map_err(|e| format!("SBDB lookup for {designation:?} failed: {e}"))?;
        let orbit = resp
            .orbit
            .ok_or_else(|| format!("SBDB returned no orbit for {designation:?}"))?;
        Ok(LiveElements {
            a: orbit.semi_major_axis,
            e: orbit.eccentricity,
            i_deg: orbit.inclination,
            omega_deg: orbit.arg_perihelion,
            omega_big_deg: orbit.long_asc_node,
            h_mag: resp.phys_par.and_then(|p| p.abs_magnitude_h),
            epoch_jd: orbit.epoch_jd,
            mean_anomaly_deg: orbit.mean_anomaly,
            mean_motion_deg_per_day: orbit.mean_motion,
            first_obs: orbit.first_obs,
            orbit_id: orbit.orbit_id,
        })
    }

    /// Diff an arbitrary snapshot table against live SBDB lookups. Errors on
    /// the first failed lookup (a missing object is a table bug, not a
    /// drift); otherwise returns every out-of-tolerance element.
    pub fn refresh_table_from_sbdb(
        client: &SbdbClient,
        snapshots: &[ElementSnapshot],
        tol: &Tolerances,
    ) -> Result<Vec<EtnoDiff>, String> {
        let mut diffs = Vec::new();
        for snapshot in snapshots {
            let live = fetch_live_elements(client, snapshot.designation)?;
            diffs.extend(diff_elements(snapshot, &live, tol));
        }
        Ok(diffs)
    }

    /// Diff `data::etno::BROWN_2017_SAMPLE` against the live SBDB with the
    /// [`Tolerances::brown2017`] tolerances. An empty result means the frozen
    /// table still matches JPL to within normal solution drift.
    pub fn refresh_etno_from_sbdb(client: &SbdbClient) -> Result<Vec<EtnoDiff>, String> {
        refresh_table_from_sbdb(client, &etno_snapshots(), &Tolerances::brown2017())
    }
}

#[cfg(feature = "sbdb-refresh")]
pub use live::{fetch_live_elements, refresh_etno_from_sbdb, refresh_table_from_sbdb};

/// Re-export so downstream crates enabling `p9-core/sbdb-refresh` can build a
/// client without depending on starfield directly.
#[cfg(feature = "sbdb-refresh")]
pub use starfield::sbdb::SbdbClient;

#[cfg(test)]
mod tests {
    use super::*;

    fn snapshot() -> ElementSnapshot {
        ElementSnapshot {
            name: "Testna",
            designation: "Testna",
            a: Some(500.0),
            e: Some(0.85),
            i_deg: Some(11.9),
            omega_deg: Some(311.5),
            omega_big_deg: Some(144.5),
            h_mag: Some(1.6),
        }
    }

    fn matching_live() -> LiveElements {
        LiveElements {
            a: Some(500.0),
            e: Some(0.85),
            i_deg: Some(11.9),
            omega_deg: Some(311.5),
            omega_big_deg: Some(144.5),
            h_mag: Some(1.6),
            ..Default::default()
        }
    }

    #[test]
    fn test_identical_elements_no_diff() {
        let diffs = diff_elements(&snapshot(), &matching_live(), &Tolerances::brown2017());
        assert!(diffs.is_empty(), "{diffs:?}");
    }

    #[test]
    fn test_drift_within_tolerance_no_diff() {
        // 0.4% in a, half the e tolerance, 0.3 deg in omega, 0.2 mag in H:
        // all within the documented Brown-2017 allowances.
        let live = LiveElements {
            a: Some(502.0),
            e: Some(0.855),
            i_deg: Some(11.9),
            omega_deg: Some(311.8),
            omega_big_deg: Some(144.5),
            h_mag: Some(1.8),
            ..Default::default()
        };
        let diffs = diff_elements(&snapshot(), &live, &Tolerances::brown2017());
        assert!(diffs.is_empty(), "{diffs:?}");
    }

    #[test]
    fn test_semi_major_axis_drift_flagged() {
        let live = LiveElements {
            a: Some(540.0), // 8% — way beyond the 1% allowance
            ..matching_live()
        };
        let diffs = diff_elements(&snapshot(), &live, &Tolerances::brown2017());
        assert_eq!(diffs.len(), 1);
        let d = &diffs[0];
        assert_eq!(d.element, Element::SemiMajorAxis);
        assert_eq!(d.object, "Testna");
        assert_eq!(d.snapshot, 500.0);
        assert_eq!(d.live, 540.0);
        assert!((d.delta - 40.0).abs() < 1e-12);
        assert!((d.tolerance - 5.0).abs() < 1e-12); // 1% of 500
    }

    #[test]
    fn test_multiple_elements_flagged_independently() {
        let live = LiveElements {
            e: Some(0.88),     // |Δe| = 0.03 > 0.01
            i_deg: Some(13.0), // 1.1 deg > 0.5
            ..matching_live()
        };
        let diffs = diff_elements(&snapshot(), &live, &Tolerances::brown2017());
        let elements: Vec<Element> = diffs.iter().map(|d| d.element).collect();
        assert_eq!(
            elements,
            vec![Element::Eccentricity, Element::Inclination],
            "{diffs:?}"
        );
    }

    #[test]
    fn test_angle_delta_is_wrap_aware() {
        // 359.9 vs 0.3 is a 0.4 deg difference, not 359.6.
        assert!((angle_delta_deg(359.9, 0.3) - 0.4).abs() < 1e-9);
        assert!((angle_delta_deg(0.3, 359.9) + 0.4).abs() < 1e-9);

        let snap = ElementSnapshot {
            omega_deg: Some(359.9),
            ..snapshot()
        };
        let live = LiveElements {
            omega_deg: Some(0.3),
            ..matching_live()
        };
        let diffs = diff_elements(&snap, &live, &Tolerances::brown2017());
        assert!(
            diffs.is_empty(),
            "wrap-around must not be flagged: {diffs:?}"
        );
    }

    #[test]
    fn test_missing_elements_skipped() {
        // A snapshot that only carries a, e, i (the neptune-crossing shape)
        // ignores live omega/Omega/H entirely.
        let snap = ElementSnapshot {
            omega_deg: None,
            omega_big_deg: None,
            h_mag: None,
            ..snapshot()
        };
        let live = LiveElements {
            omega_deg: Some(100.0),    // wildly different but not compared
            omega_big_deg: Some(10.0), // ditto
            h_mag: None,
            ..matching_live()
        };
        let diffs = diff_elements(&snap, &live, &Tolerances::brown2017());
        assert!(diffs.is_empty(), "{diffs:?}");
    }

    #[test]
    fn test_etno_snapshots_cover_table() {
        let snaps = etno_snapshots();
        assert_eq!(snaps.len(), BROWN_2017_SAMPLE.len());
        assert!(snaps.iter().all(|s| s.a.is_some()
            && s.e.is_some()
            && s.i_deg.is_some()
            && s.omega_deg.is_some()
            && s.omega_big_deg.is_some()
            && s.h_mag.is_some()));
        // A snapshot diffed against itself is clean.
        for s in &snaps {
            let live = LiveElements {
                a: s.a,
                e: s.e,
                i_deg: s.i_deg,
                omega_deg: s.omega_deg,
                omega_big_deg: s.omega_big_deg,
                h_mag: s.h_mag,
                ..Default::default()
            };
            assert!(diff_elements(s, &live, &Tolerances::brown2017()).is_empty());
        }
    }

    #[test]
    fn test_calendar_to_jd_known_epochs() {
        // J2000.0
        assert!((calendar_to_jd(2000, 1, 1.5) - 2_451_545.0).abs() < 1e-9);
        // 1987-04-10.0 = JD 2446895.5 (Meeus, example 7.b)
        assert!((calendar_to_jd(1987, 4, 10.0) - 2_446_895.5).abs() < 1e-9);
    }

    #[test]
    fn test_parse_sbdb_date() {
        assert!((parse_sbdb_date("2000-01-01").unwrap() - 2_451_544.5).abs() < 1e-9);
        assert!((parse_sbdb_date("2013-04-04").unwrap() - 2_456_386.5).abs() < 1e-9);
        assert_eq!(parse_sbdb_date("garbage"), None);
        assert_eq!(parse_sbdb_date("2013-13-04"), None);
    }

    #[test]
    fn test_heliocentric_distance_at_epoch_and_extremes() {
        // a = 100 AU, e = 0.5, at epoch M0 = 0 (perihelion): r = q = 50.
        let live = LiveElements {
            a: Some(100.0),
            e: Some(0.5),
            epoch_jd: Some(2_451_545.0),
            mean_anomaly_deg: Some(0.0),
            mean_motion_deg_per_day: Some(360.0 / 365_250.0), // 1000 yr period
            ..Default::default()
        };
        let r_peri = heliocentric_distance_at(&live, 2_451_545.0).unwrap();
        assert!((r_peri - 50.0).abs() < 1e-9, "r = {r_peri}");
        // Half a period later: aphelion, r = a(1+e) = 150.
        let r_apo = heliocentric_distance_at(&live, 2_451_545.0 + 365_250.0 / 2.0).unwrap();
        assert!((r_apo - 150.0).abs() < 1e-6, "r = {r_apo}");
        // Always within [q, Q].
        for k in 0..20 {
            let jd = 2_451_545.0 + k as f64 * 20_000.0;
            let r = heliocentric_distance_at(&live, jd).unwrap();
            assert!((50.0..=150.0).contains(&r), "r = {r}");
        }
    }

    #[test]
    fn test_first_obs_distance_requires_fields() {
        let mut live = LiveElements {
            a: Some(100.0),
            e: Some(0.5),
            epoch_jd: Some(2_451_545.0),
            mean_anomaly_deg: Some(0.0),
            mean_motion_deg_per_day: Some(0.001),
            first_obs: None,
            ..Default::default()
        };
        assert_eq!(first_obs_distance(&live), None);
        live.first_obs = Some("2000-01-01".into());
        let r = first_obs_distance(&live).unwrap();
        assert!((50.0..=150.0).contains(&r));
    }

    #[test]
    fn test_diff_display_names_object_and_element() {
        let live = LiveElements {
            a: Some(540.0),
            ..matching_live()
        };
        let diffs = diff_elements(&snapshot(), &live, &Tolerances::brown2017());
        let msg = diffs[0].to_string();
        assert!(msg.contains("Testna"), "{msg}");
        assert!(msg.contains("a [AU]"), "{msg}");
        assert!(msg.contains("540"), "{msg}");
    }
}
