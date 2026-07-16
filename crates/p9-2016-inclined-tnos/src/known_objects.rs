//! Known highly inclined TNOs from the paper.
//!
//! These are the observed objects that the simulation aims to reproduce:
//! - Drac (2008 KV42): retrograde centaur, i ≈ 103°
//! - Niku (2011 KT19): high-i TNO, i ≈ 110°
//! - 2016 NM56: retrograde, i ≈ 144°

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;

/// A known high-inclination TNO for comparison.
#[derive(Debug, Clone)]
pub struct KnownTno {
    pub name: &'static str,
    pub designation: &'static str,
    pub elements: OrbitalElements,
}

/// Return the three key high-inclination TNOs discussed in the paper.
///
/// Elements re-transcribed from JPL SBDB 2026-07. The previous table carried
/// wrong angular elements for four of five objects (e.g. Drac's ω duplicated
/// its node; Niku's ω off by 240°) and a "2013 LA2" row matching no real
/// object — the real 2013 LA2 is a ~5.7 AU retrograde centaur, so that row
/// was dropped.
pub fn paper_tnos() -> Vec<KnownTno> {
    vec![
        KnownTno {
            name: "Drac",
            designation: "2008 KV42",
            elements: OrbitalElements {
                a: 41.5,
                e: 0.491,
                i: 103.0 * DEG2RAD,
                omega: 133.0 * DEG2RAD,
                omega_big: 261.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
        KnownTno {
            name: "Niku",
            designation: "2011 KT19",
            elements: OrbitalElements {
                a: 35.6,
                e: 0.333,
                i: 110.0 * DEG2RAD,
                omega: 322.0 * DEG2RAD,
                omega_big: 244.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
        KnownTno {
            name: "2016 NM56",
            designation: "2016 NM56",
            elements: OrbitalElements {
                a: 73.3,
                e: 0.856,
                i: 144.0 * DEG2RAD,
                omega: 346.0 * DEG2RAD,
                omega_big: 350.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
        },
    ]
}

/// Additional high-inclination centaurs and TNOs for extended comparison
/// (JPL SBDB 2026-07 elements).
pub fn extended_high_i_objects() -> Vec<KnownTno> {
    let mut objects = paper_tnos();
    objects.push(KnownTno {
        name: "2010 WG9",
        designation: "2010 WG9",
        elements: OrbitalElements {
            a: 53.5,
            e: 0.648,
            i: 70.2 * DEG2RAD,
            omega: 293.0 * DEG2RAD,
            omega_big: 92.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
    });
    objects
}

/// The paper table as snapshot rows for `p9_core::data::refresh` (all five
/// elements are carried).
pub fn element_snapshots() -> Vec<p9_core::data::refresh::ElementSnapshot> {
    paper_tnos()
        .iter()
        .map(|t| p9_core::data::refresh::ElementSnapshot {
            name: t.name,
            designation: t.designation,
            a: Some(t.elements.a),
            e: Some(t.elements.e),
            i_deg: Some(t.elements.i / DEG2RAD),
            omega_deg: Some(t.elements.omega / DEG2RAD),
            omega_big_deg: Some(t.elements.omega_big / DEG2RAD),
            h_mag: None,
        })
        .collect()
}

/// Diff the frozen table against the live JPL SBDB (network; see the
/// `#[ignore]`d test in `tests/sbdb_refresh_live.rs`). The table was
/// re-transcribed at full SBDB precision in July 2026, so the tight
/// `full_precision` tolerances apply.
#[cfg(feature = "sbdb-refresh")]
pub fn refresh_from_sbdb(
    client: &p9_core::data::refresh::SbdbClient,
) -> Result<Vec<p9_core::data::refresh::EtnoDiff>, String> {
    p9_core::data::refresh::refresh_table_from_sbdb(
        client,
        &element_snapshots(),
        &p9_core::data::refresh::Tolerances::full_precision(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paper_tnos() {
        let tnos = paper_tnos();
        assert_eq!(tnos.len(), 3);

        for tno in &tnos {
            // All should have i > 90° (retrograde or near-retrograde)
            assert!(
                tno.elements.i > 90.0 * DEG2RAD,
                "{} has i = {:.1}° (expected > 90°)",
                tno.name,
                tno.elements.i / DEG2RAD
            );
            // All should have a < 100 AU (decoupled from P9)
            assert!(
                tno.elements.a < 100.0,
                "{} has a = {:.1} AU (expected < 100)",
                tno.name,
                tno.elements.a
            );
        }
    }

    #[test]
    fn test_drac_is_retrograde() {
        let tnos = paper_tnos();
        let drac = &tnos[0];
        assert!(drac.elements.i > 90.0 * DEG2RAD);
        assert!((drac.elements.a - 41.4).abs() < 1.0);
    }
}
