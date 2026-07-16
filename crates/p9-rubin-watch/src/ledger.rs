//! The distant-object ledger: the committed record of every object the watch
//! tracks, with provenance, quality, class tags, and an append-only history.
//!
//! Schema v1 per `rubin_watch/design/07-schemas.md`. Serialization is
//! canonical (fixed struct field order, objects sorted by id) and the file is
//! rewritten only when content changes, so a no-op sync leaves the ledger
//! byte-identical (the Phase-1 idempotence requirement).

use serde::{Deserialize, Serialize};

pub const SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Ledger {
    pub schema_version: u32,
    pub generated_utc: String,
    pub source_stamps: SourceStamps,
    pub objects: Vec<ObjectRecord>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SourceStamps {
    pub mpc_distant_last_modified: Option<String>,
    pub sbdb_query_fetched: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ObjectRecord {
    /// Stable internal key, "rw-NNNN"; never reused, never re-keyed.
    pub id: String,
    pub primary_desig: String,
    pub aliases: Vec<String>,
    pub tier: Tier,
    pub class_tags: Vec<String>,
    pub elements: Elements,
    /// Populated when MPC and SBDB disagree beyond the refresh tolerances.
    pub elements_alt: Option<Elements>,
    pub quality: Quality,
    pub flags: Vec<String>,
    pub history: Vec<HistoryEvent>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Tier {
    Provisional,
    Vetted,
    Retired,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Elements {
    pub epoch_jd: f64,
    pub a_au: f64,
    pub e: f64,
    pub q_au: f64,
    pub i_deg: f64,
    pub omega_deg: f64,
    pub omega_big_deg: f64,
    pub varpi_deg: f64,
    pub h_mag: Option<f64>,
    /// "sbdb" | "mpc_distant"
    pub source: String,
    pub fetched_utc: String,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Quality {
    pub arc_days: Option<f64>,
    pub n_oppositions: Option<u32>,
    /// MPC U parameter / SBDB condition code, 0 (best) .. 9; None = unknown.
    pub condition_code: Option<u32>,
    pub n_obs: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HistoryEvent {
    pub utc: String,
    /// ingested | refreshed | promoted | demoted | reclassified | renamed | retired
    pub event: String,
    pub detail: String,
}

impl Ledger {
    pub fn empty(now_utc: &str) -> Self {
        Ledger {
            schema_version: SCHEMA_VERSION,
            generated_utc: now_utc.to_string(),
            source_stamps: SourceStamps::default(),
            objects: Vec::new(),
        }
    }

    pub fn load(path: &std::path::Path) -> Result<Self, String> {
        let raw = std::fs::read_to_string(path).map_err(|e| format!("read {path:?}: {e}"))?;
        let ledger: Ledger =
            serde_json::from_str(&raw).map_err(|e| format!("parse {path:?}: {e}"))?;
        if ledger.schema_version != SCHEMA_VERSION {
            return Err(format!(
                "ledger schema_version {} != supported {SCHEMA_VERSION}",
                ledger.schema_version
            ));
        }
        ledger.check_invariants()?;
        Ok(ledger)
    }

    /// Canonical serialization: objects sorted by id, pretty JSON, trailing
    /// newline.
    pub fn to_canonical_json(&self) -> String {
        let mut sorted = self.clone();
        sorted.objects.sort_by(|a, b| a.id.cmp(&b.id));
        let mut s = serde_json::to_string_pretty(&sorted).expect("ledger serializes");
        s.push('\n');
        s
    }

    /// Write only if the content (ignoring `generated_utc`) differs from what
    /// is on disk. Returns true when a write happened.
    pub fn save_if_changed(&self, path: &std::path::Path) -> Result<bool, String> {
        if let Ok(existing) = Ledger::load(path) {
            if content_equal(&existing, self) {
                return Ok(false);
            }
        }
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir).map_err(|e| format!("mkdir {dir:?}: {e}"))?;
        }
        std::fs::write(path, self.to_canonical_json()).map_err(|e| format!("write: {e}"))?;
        Ok(true)
    }

    /// Resolve any known designation (primary or alias) to the object index.
    pub fn find_by_desig(&self, desig: &str) -> Option<usize> {
        self.objects
            .iter()
            .position(|o| o.primary_desig == desig || o.aliases.iter().any(|a| a == desig))
    }

    pub fn next_id(&self) -> String {
        let max: u32 = self
            .objects
            .iter()
            .filter_map(|o| o.id.strip_prefix("rw-").and_then(|n| n.parse().ok()))
            .max()
            .unwrap_or(0);
        format!("rw-{:04}", max + 1)
    }

    pub fn check_invariants(&self) -> Result<(), String> {
        let mut ids = std::collections::HashSet::new();
        let mut aliases = std::collections::HashMap::new();
        for o in &self.objects {
            if !ids.insert(&o.id) {
                return Err(format!("duplicate id {}", o.id));
            }
            for a in std::iter::once(&o.primary_desig).chain(o.aliases.iter()) {
                if let Some(prev) = aliases.insert(a.clone(), o.id.clone()) {
                    if prev != o.id {
                        return Err(format!("alias {a:?} maps to both {prev} and {}", o.id));
                    }
                }
            }
            if o.history.is_empty() {
                return Err(format!("{}: empty history", o.id));
            }
        }
        Ok(())
    }
}

/// Equality ignoring the volatile stamps (`generated_utc`, `source_stamps`):
/// the ledger file rewrites only when the OBJECTS change.
pub fn content_equal(a: &Ledger, b: &Ledger) -> bool {
    let mut a2 = a.clone();
    let mut b2 = b.clone();
    a2.generated_utc = String::new();
    b2.generated_utc = String::new();
    a2.source_stamps = SourceStamps::default();
    b2.source_stamps = SourceStamps::default();
    a2.objects.sort_by(|x, y| x.id.cmp(&y.id));
    b2.objects.sort_by(|x, y| x.id.cmp(&y.id));
    a2 == b2
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(id: &str, desig: &str) -> ObjectRecord {
        ObjectRecord {
            id: id.into(),
            primary_desig: desig.into(),
            aliases: vec![desig.into()],
            tier: Tier::Provisional,
            class_tags: vec!["watch_only".into()],
            elements: Elements {
                epoch_jd: 2461200.5,
                a_au: 200.0,
                e: 0.5,
                q_au: 100.0,
                i_deg: 10.0,
                omega_deg: 0.0,
                omega_big_deg: 0.0,
                varpi_deg: 0.0,
                h_mag: None,
                source: "sbdb".into(),
                fetched_utc: "2026-07-16T00:00:00Z".into(),
            },
            elements_alt: None,
            quality: Quality::default(),
            flags: vec![],
            history: vec![HistoryEvent {
                utc: "2026-07-16T00:00:00Z".into(),
                event: "ingested".into(),
                detail: "test".into(),
            }],
        }
    }

    #[test]
    fn canonical_serialization_is_stable_and_sorted() {
        let mut l = Ledger::empty("2026-07-16T00:00:00Z");
        l.objects.push(record("rw-0002", "B"));
        l.objects.push(record("rw-0001", "A"));
        let s1 = l.to_canonical_json();
        let s2 = l.to_canonical_json();
        assert_eq!(s1, s2);
        assert!(s1.find("rw-0001").unwrap() < s1.find("rw-0002").unwrap());
        assert!(s1.ends_with('\n'));
    }

    #[test]
    fn save_if_changed_is_idempotent() {
        let dir = std::env::temp_dir().join("p9rw-ledger-test");
        let path = dir.join("ledger.json");
        let _ = std::fs::remove_file(&path);
        let mut l = Ledger::empty("2026-07-16T00:00:00Z");
        l.objects.push(record("rw-0001", "A"));
        assert!(l.save_if_changed(&path).unwrap());
        let bytes1 = std::fs::read(&path).unwrap();
        // Same content, different stamp: must NOT rewrite.
        let mut l2 = l.clone();
        l2.generated_utc = "2026-07-17T00:00:00Z".into();
        assert!(!l2.save_if_changed(&path).unwrap());
        assert_eq!(bytes1, std::fs::read(&path).unwrap());
        // Real change: rewrites.
        l2.objects.push(record("rw-0002", "B"));
        assert!(l2.save_if_changed(&path).unwrap());
    }

    #[test]
    fn invariants_catch_duplicate_alias() {
        let mut l = Ledger::empty("2026-07-16T00:00:00Z");
        l.objects.push(record("rw-0001", "A"));
        let mut b = record("rw-0002", "B");
        b.aliases.push("A".into());
        l.objects.push(b);
        assert!(l.check_invariants().is_err());
    }

    #[test]
    fn find_resolves_aliases_and_next_id_advances() {
        let mut l = Ledger::empty("2026-07-16T00:00:00Z");
        let mut r = record("rw-0007", "2003 VB12");
        r.aliases.push("Sedna".into());
        l.objects.push(r);
        assert_eq!(l.find_by_desig("Sedna"), Some(0));
        assert_eq!(l.find_by_desig("2003 VB12"), Some(0));
        assert_eq!(l.find_by_desig("nope"), None);
        assert_eq!(l.next_id(), "rw-0008");
    }
}
