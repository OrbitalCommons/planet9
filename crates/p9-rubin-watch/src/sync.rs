//! The pure sync step: (ledger, merged census, now) → (updated ledger, delta).
//!
//! No I/O here — everything is fixture-testable. Semantics per design/02:
//! new objects are ingested; known objects are refreshed, re-classified and
//! re-tiered with history events; objects whose revised orbits fall out of
//! the census cut are retired (never deleted).

use crate::classify::{class_tags, in_census, tier};
use crate::ledger::{Elements, HistoryEvent, Ledger, ObjectRecord, Tier};
use crate::sources::MergedRecord;

#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct Delta {
    pub new: Vec<String>,
    pub promoted: Vec<String>,
    pub demoted: Vec<String>,
    pub reclassified: Vec<String>,
    pub renamed: Vec<String>,
    pub retired: Vec<String>,
    pub discrepant: Vec<String>,
    pub refreshed: usize,
}

impl Delta {
    pub fn is_empty(&self) -> bool {
        self.new.is_empty()
            && self.promoted.is_empty()
            && self.demoted.is_empty()
            && self.reclassified.is_empty()
            && self.renamed.is_empty()
            && self.retired.is_empty()
            && self.discrepant.is_empty()
    }
}

fn to_elements(r: &crate::sources::CensusRecord, now_utc: &str) -> Elements {
    Elements {
        epoch_jd: r.epoch_jd,
        a_au: r.a_au,
        e: r.e,
        q_au: r.q_au,
        i_deg: r.i_deg,
        omega_deg: r.omega_deg,
        omega_big_deg: r.omega_big_deg,
        varpi_deg: r.varpi_deg(),
        h_mag: r.h_mag,
        source: r.source.clone(),
        fetched_utc: now_utc.to_string(),
    }
}

/// Element change big enough to record a refresh (below this, leave the
/// ledger untouched so no-op syncs stay byte-identical).
fn materially_different(old: &Elements, new: &Elements) -> bool {
    (old.a_au - new.a_au).abs() > 1e-9
        || (old.e - new.e).abs() > 1e-12
        || (old.q_au - new.q_au).abs() > 1e-9
        || (old.i_deg - new.i_deg).abs() > 1e-9
        || (old.omega_deg - new.omega_deg).abs() > 1e-9
        || (old.omega_big_deg - new.omega_big_deg).abs() > 1e-9
}

pub fn sync(mut ledger: Ledger, census: &[MergedRecord], now_utc: &str) -> (Ledger, Delta) {
    let mut delta = Delta::default();

    for m in census {
        let r = &m.best;
        let known = r.aliases.iter().find_map(|a| ledger.find_by_desig(a));
        match known {
            None => {
                let id = ledger.next_id();
                let quality = crate::ledger::Quality {
                    arc_days: r.arc_days,
                    n_oppositions: r.n_oppositions,
                    condition_code: r.condition_code,
                    n_obs: r.n_obs,
                };
                let mut flags = Vec::new();
                if m.alt.is_some() {
                    flags.push("discrepant".to_string());
                    delta.discrepant.push(r.desig.clone());
                }
                ledger.objects.push(ObjectRecord {
                    id: id.clone(),
                    primary_desig: r.desig.clone(),
                    aliases: r.aliases.clone(),
                    tier: tier(&quality, r.e),
                    class_tags: class_tags(r.a_au, r.q_au, r.i_deg),
                    elements: to_elements(r, now_utc),
                    elements_alt: m.alt.as_ref().map(|a| to_elements(a, now_utc)),
                    quality,
                    flags,
                    history: vec![HistoryEvent {
                        utc: now_utc.to_string(),
                        event: "ingested".into(),
                        detail: format!("census {}", r.source),
                    }],
                });
                delta.new.push(r.desig.clone());
            }
            Some(idx) => {
                let obj = &mut ledger.objects[idx];
                // New aliases (numbering / naming events).
                let mut renamed = false;
                for a in &r.aliases {
                    if !obj.aliases.contains(a) {
                        obj.aliases.push(a.clone());
                        renamed = true;
                    }
                }
                if renamed {
                    obj.history.push(HistoryEvent {
                        utc: now_utc.to_string(),
                        event: "renamed".into(),
                        detail: format!("aliases now {:?}", obj.aliases),
                    });
                    delta.renamed.push(obj.primary_desig.clone());
                }

                // Element refresh.
                let new_el = to_elements(r, now_utc);
                if materially_different(&obj.elements, &new_el) {
                    obj.elements = new_el;
                    obj.quality = crate::ledger::Quality {
                        arc_days: r.arc_days,
                        n_oppositions: r.n_oppositions,
                        condition_code: r.condition_code,
                        n_obs: r.n_obs,
                    };
                    obj.history.push(HistoryEvent {
                        utc: now_utc.to_string(),
                        event: "refreshed".into(),
                        detail: format!("elements from {}", r.source),
                    });
                    delta.refreshed += 1;
                }

                // Discrepancy flag tracks the current merge.
                let flagged = obj.flags.iter().any(|f| f == "discrepant");
                match (&m.alt, flagged) {
                    (Some(alt), _) => {
                        let new_alt = to_elements(alt, now_utc);
                        let alt_changed = obj
                            .elements_alt
                            .as_ref()
                            .is_none_or(|old| materially_different(old, &new_alt));
                        if alt_changed {
                            obj.elements_alt = Some(new_alt);
                        }
                        if !flagged {
                            obj.flags.push("discrepant".into());
                            delta.discrepant.push(obj.primary_desig.clone());
                        }
                    }
                    (None, true) => {
                        obj.flags.retain(|f| f != "discrepant");
                        obj.elements_alt = None;
                    }
                    (None, false) => {}
                }

                // Re-classify and re-tier.
                let tags = class_tags(r.a_au, r.q_au, r.i_deg);
                if tags != obj.class_tags {
                    obj.history.push(HistoryEvent {
                        utc: now_utc.to_string(),
                        event: "reclassified".into(),
                        detail: format!("{:?} -> {:?}", obj.class_tags, tags),
                    });
                    obj.class_tags = tags;
                    delta.reclassified.push(obj.primary_desig.clone());
                }
                let new_tier = tier(&obj.quality, r.e);
                if new_tier != obj.tier {
                    let (event, list) = if new_tier == Tier::Vetted {
                        ("promoted", &mut delta.promoted)
                    } else {
                        ("demoted", &mut delta.demoted)
                    };
                    obj.history.push(HistoryEvent {
                        utc: now_utc.to_string(),
                        event: event.into(),
                        detail: format!("{:?} -> {:?}", obj.tier, new_tier),
                    });
                    obj.tier = new_tier;
                    list.push(obj.primary_desig.clone());
                }
            }
        }
    }

    // Retirement: a tracked object whose *revised* orbit no longer passes the
    // census cut (element revisions arrive through the loop above first).
    for obj in ledger.objects.iter_mut() {
        if obj.tier != Tier::Retired
            && !in_census(obj.elements.a_au, obj.elements.q_au, obj.elements.i_deg)
        {
            obj.history.push(HistoryEvent {
                utc: now_utc.to_string(),
                event: "retired".into(),
                detail: "orbit revision fell below census cut".into(),
            });
            obj.tier = Tier::Retired;
            delta.retired.push(obj.primary_desig.clone());
        }
    }

    ledger.generated_utc = now_utc.to_string();
    (ledger, delta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources::CensusRecord;

    fn rec(desig: &str, a: f64, e: f64, q: f64, i: f64, n_opp: u32, cc: u32) -> MergedRecord {
        MergedRecord {
            best: CensusRecord {
                desig: desig.into(),
                aliases: vec![desig.into()],
                epoch_jd: 2461200.5,
                a_au: a,
                e,
                q_au: q,
                i_deg: i,
                omega_deg: 100.0,
                omega_big_deg: 50.0,
                h_mag: Some(5.0),
                arc_days: Some(5000.0),
                n_oppositions: Some(n_opp),
                condition_code: Some(cc),
                n_obs: Some(100),
                source: "sbdb".into(),
            },
            alt: None,
        }
    }

    #[test]
    fn new_object_ingests_classifies_and_tiers() {
        let ledger = Ledger::empty("t0");
        let census = vec![rec("2099 XX1", 400.0, 0.88, 48.0, 15.0, 1, 9)];
        let (ledger, delta) = sync(ledger, &census, "t1");
        assert_eq!(delta.new, vec!["2099 XX1"]);
        let o = &ledger.objects[0];
        assert_eq!(o.id, "rw-0001");
        assert_eq!(o.tier, Tier::Provisional); // 1 opp, U=9
        assert!(o.class_tags.contains(&"etno_vetted_class".into()));
        assert_eq!(o.elements.varpi_deg, 150.0);
    }

    #[test]
    fn idempotent_second_sync_is_a_noop() {
        let census = vec![rec("2099 XX1", 400.0, 0.88, 48.0, 15.0, 3, 2)];
        let (ledger, _) = sync(Ledger::empty("t0"), &census, "t1");
        let before = ledger.to_canonical_json();
        let (ledger2, delta) = sync(ledger, &census, "t2");
        assert!(delta.is_empty(), "{delta:?}");
        // generated_utc moves but content is identical.
        assert!(crate::ledger::content_equal(
            &serde_json::from_str(&before).unwrap(),
            &ledger2
        ));
    }

    #[test]
    fn rename_appends_alias_without_duplicate() {
        let census = vec![rec("2099 XX1", 400.0, 0.88, 48.0, 15.0, 3, 2)];
        let (ledger, _) = sync(Ledger::empty("t0"), &census, "t1");
        let mut renamed = rec("2099 XX1", 400.0, 0.88, 48.0, 15.0, 3, 2);
        renamed.best.aliases.push("(612345)".into());
        let (ledger, delta) = sync(ledger, &[renamed], "t2");
        assert_eq!(ledger.objects.len(), 1);
        assert_eq!(delta.renamed.len(), 1);
        assert!(ledger.objects[0].aliases.contains(&"(612345)".into()));
    }

    #[test]
    fn promotion_and_demotion_paths() {
        // Provisional first (1 opp)...
        let (ledger, _) = sync(
            Ledger::empty("t0"),
            &[rec("2099 XX1", 400.0, 0.88, 48.0, 15.0, 1, 9)],
            "t1",
        );
        assert_eq!(ledger.objects[0].tier, Tier::Provisional);
        // ...second opposition with a good fit promotes...
        let (ledger, delta) = sync(
            ledger,
            &[rec("2099 XX1", 405.0, 0.881, 48.2, 15.0, 2, 3)],
            "t2",
        );
        assert_eq!(delta.promoted.len(), 1);
        assert_eq!(ledger.objects[0].tier, Tier::Vetted);
        // ...a shocking re-fit demotes.
        let (ledger, delta) = sync(
            ledger,
            &[rec("2099 XX1", 300.0, 0.84, 48.0, 15.0, 2, 8)],
            "t3",
        );
        assert_eq!(delta.demoted.len(), 1);
        assert_eq!(ledger.objects[0].tier, Tier::Provisional);
    }

    #[test]
    fn orbit_revision_below_cut_retires_and_reclassifies() {
        let (ledger, _) = sync(
            Ledger::empty("t0"),
            &[rec("2099 XX1", 260.0, 0.83, 44.0, 15.0, 2, 2)],
            "t1",
        );
        assert!(ledger.objects[0]
            .class_tags
            .contains(&"etno_vetted_class".into()));
        // Revision drops it to a classical-belt orbit: retired.
        let (ledger, delta) = sync(
            ledger,
            &[rec("2099 XX1", 45.0, 0.05, 42.7, 5.0, 2, 2)],
            "t2",
        );
        assert_eq!(delta.retired.len(), 1);
        assert_eq!(ledger.objects[0].tier, Tier::Retired);
        assert_eq!(delta.reclassified.len(), 1);
    }
}
