//! Census sources: the MPC distant-object extended file and the JPL SBDB
//! query API, parsed into a common [`CensusRecord`] and merged.
//!
//! Parsers are pure (fixture-tested offline); the two `fetch_*` functions at
//! the bottom are the only network touchers in the crate and run only from
//! the `mpc-sync` subcommand (SR-8).

use crate::classify::in_census;
use serde_json::Value;

pub const MPC_DISTANT_URL: &str =
    "https://www.minorplanetcenter.net/Extended_Files/distant_extended.json.gz";
pub const SBDB_QUERY_URL: &str = "https://ssd-api.jpl.nasa.gov/sbdb_query.api";
/// Minor planets only (`sb-kind=a`): long-period comets flood the a > 100 AU
/// cut and their discoverer "names" (Pons, Tempel, "Great comet") collide
/// across objects, poisoning the alias join. The MPC distant file is
/// minor-planets-only, so this also keeps the two sources congruent.
pub const SBDB_KIND: &str = "a";
/// The census constraint, mirroring `classify::in_census` server-side.
pub const SBDB_CDATA: &str =
    r#"{"OR":[{"AND":["a|GT|100"]},{"AND":["i|GT|60","a|GT|30","q|GT|15"]}]}"#;
pub const SBDB_FIELDS: &str =
    "pdes,name,full_name,a,e,q,i,om,w,epoch,condition_code,data_arc,n_opp,n_obs_used,H";

/// One object as seen by one source, normalized.
#[derive(Debug, Clone, PartialEq)]
pub struct CensusRecord {
    /// Join key: provisional designation when known (e.g. "2003 VB12"),
    /// else number or name.
    pub desig: String,
    /// Every designation string seen for this object (includes `desig`).
    pub aliases: Vec<String>,
    pub epoch_jd: f64,
    pub a_au: f64,
    pub e: f64,
    pub q_au: f64,
    pub i_deg: f64,
    pub omega_deg: f64,
    pub omega_big_deg: f64,
    pub h_mag: Option<f64>,
    pub arc_days: Option<f64>,
    pub n_oppositions: Option<u32>,
    pub condition_code: Option<u32>,
    pub n_obs: Option<u32>,
    /// "sbdb" | "mpc_distant"
    pub source: String,
}

impl CensusRecord {
    pub fn varpi_deg(&self) -> f64 {
        (self.omega_deg + self.omega_big_deg).rem_euclid(360.0)
    }
}

// ---------------------------------------------------------------- MPC parse

/// Parse the distant-object extended file (already gunzipped), keeping only
/// census members.
pub fn parse_mpc_distant(json: &str) -> Result<Vec<CensusRecord>, String> {
    let rows: Vec<Value> = serde_json::from_str(json).map_err(|e| format!("mpc json: {e}"))?;
    let mut out = Vec::new();
    for r in &rows {
        let f = |k: &str| r.get(k).and_then(Value::as_f64);
        let (Some(a), Some(e), Some(i)) = (f("a"), f("e"), f("i")) else {
            return Err(format!("mpc row missing a/e/i: {r}"));
        };
        let q = f("Perihelion_dist").unwrap_or(a * (1.0 - e));
        if !in_census(a, q, i) {
            continue;
        }
        let s = |k: &str| r.get(k).and_then(Value::as_str).map(str::to_string);
        let principal = s("Principal_desig");
        let name = s("Name");
        let number = s("Number").map(|n| n.trim_matches(['(', ')']).to_string());
        let desig = principal
            .clone()
            .or_else(|| number.clone())
            .or_else(|| name.clone())
            .ok_or_else(|| format!("mpc row with no designation: {r}"))?;
        let mut aliases = vec![desig.clone()];
        for extra in [principal, name, number].into_iter().flatten() {
            if !aliases.contains(&extra) {
                aliases.push(extra);
            }
        }
        out.push(CensusRecord {
            desig,
            aliases,
            epoch_jd: f("Epoch").unwrap_or(0.0),
            a_au: a,
            e,
            q_au: q,
            i_deg: i,
            omega_deg: f("Peri").unwrap_or(0.0),
            omega_big_deg: f("Node").unwrap_or(0.0),
            h_mag: f("H"),
            arc_days: mpc_arc_days(r),
            n_oppositions: f("Num_opps").map(|n| n as u32),
            condition_code: r
                .get("U")
                .and_then(Value::as_str)
                .and_then(|u| u.trim().parse::<u32>().ok()),
            n_obs: f("Num_obs").map(|n| n as u32),
            source: "mpc_distant".into(),
        });
    }
    Ok(out)
}

/// Arc length in days: `Arc_length` (single-opposition, days) or the span of
/// `Arc_years` ("1990-2026").
fn mpc_arc_days(r: &Value) -> Option<f64> {
    if let Some(d) = r.get("Arc_length").and_then(Value::as_f64) {
        return Some(d);
    }
    let years = r.get("Arc_years").and_then(Value::as_str)?;
    let (y1, y2) = years.split_once('-')?;
    let (y1, y2): (f64, f64) = (y1.trim().parse().ok()?, y2.trim().parse().ok()?);
    Some((y2 - y1).max(0.0) * 365.25)
}

// --------------------------------------------------------------- SBDB parse

/// Parse an sbdb_query.api response (fields per [`SBDB_FIELDS`]).
pub fn parse_sbdb_query(json: &str) -> Result<Vec<CensusRecord>, String> {
    let v: Value = serde_json::from_str(json).map_err(|e| format!("sbdb json: {e}"))?;
    let fields: Vec<String> = v
        .get("fields")
        .and_then(Value::as_array)
        .ok_or("sbdb: no fields")?
        .iter()
        .filter_map(|f| f.as_str().map(str::to_string))
        .collect();
    let idx = |name: &str| fields.iter().position(|f| f == name);
    let (ia, ie, iq, ii) = (
        idx("a").ok_or("no a")?,
        idx("e").ok_or("no e")?,
        idx("q").ok_or("no q")?,
        idx("i").ok_or("no i")?,
    );
    let data = v
        .get("data")
        .and_then(Value::as_array)
        .ok_or("sbdb: no data")?;
    let cell_f = |row: &[Value], i: usize| -> Option<f64> {
        row.get(i).and_then(|c| match c {
            Value::Number(n) => n.as_f64(),
            Value::String(s) => s.trim().parse().ok(),
            _ => None,
        })
    };
    let cell_s = |row: &[Value], i: Option<usize>| -> Option<String> {
        i.and_then(|i| row.get(i))
            .and_then(Value::as_str)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
    };
    let mut out = Vec::new();
    for row in data {
        let row = row.as_array().ok_or("sbdb row not array")?;
        let (Some(a), Some(e), Some(q), Some(i)) = (
            cell_f(row, ia),
            cell_f(row, ie),
            cell_f(row, iq),
            cell_f(row, ii),
        ) else {
            continue; // unconstrained orbit rows can miss elements
        };
        if !in_census(a, q, i) {
            continue;
        }
        let pdes = cell_s(row, idx("pdes"));
        let name = cell_s(row, idx("name"));
        let full = cell_s(row, idx("full_name"));
        // Provisional designation lives in full_name's parenthetical for
        // numbered objects: "90377 Sedna (2003 VB12)".
        let provisional = full
            .as_deref()
            .and_then(|f| {
                f.rsplit_once('(')
                    .map(|(_, t)| t.trim_end_matches(')').trim())
            })
            .filter(|p| looks_provisional(p))
            .map(str::to_string);
        let desig = provisional
            .clone()
            .or_else(|| pdes.clone())
            .or_else(|| name.clone())
            .ok_or("sbdb row with no designation")?;
        // A bare name is only a safe alias for NUMBERED minor planets (IAU
        // names are unique there); unnumbered/cometary name fields are not.
        let numbered = pdes
            .as_deref()
            .is_some_and(|p| p.chars().all(|c| c.is_ascii_digit()));
        let name_alias = if numbered { name.clone() } else { None };
        let mut aliases = vec![desig.clone()];
        for extra in [provisional, pdes, name_alias].into_iter().flatten() {
            if !aliases.contains(&extra) {
                aliases.push(extra);
            }
        }
        out.push(CensusRecord {
            desig,
            aliases,
            epoch_jd: cell_f(row, idx("epoch").ok_or("no epoch")?).unwrap_or(0.0),
            a_au: a,
            e,
            q_au: q,
            i_deg: i,
            omega_deg: idx("w").and_then(|i| cell_f(row, i)).unwrap_or(0.0),
            omega_big_deg: idx("om").and_then(|i| cell_f(row, i)).unwrap_or(0.0),
            h_mag: idx("H").and_then(|i| cell_f(row, i)),
            arc_days: idx("data_arc").and_then(|i| cell_f(row, i)),
            n_oppositions: idx("n_opp").and_then(|i| cell_f(row, i)).map(|n| n as u32),
            condition_code: idx("condition_code")
                .and_then(|i| cell_f(row, i))
                .map(|c| c as u32),
            n_obs: idx("n_obs_used")
                .and_then(|i| cell_f(row, i))
                .map(|n| n as u32),
            source: "sbdb".into(),
        });
    }
    Ok(out)
}

/// "2003 VB12"-shaped?
fn looks_provisional(s: &str) -> bool {
    let mut parts = s.split_whitespace();
    let (Some(y), Some(code), None) = (parts.next(), parts.next(), parts.next()) else {
        return false;
    };
    y.len() == 4
        && y.chars().all(|c| c.is_ascii_digit())
        && code.chars().next().is_some_and(|c| c.is_ascii_uppercase())
}

// -------------------------------------------------------------------- merge

/// Merge the two source lists into one census keyed by designation: SBDB
/// elements win (full precision), MPC supplies opposition counts / U where
/// SBDB lacks them, and material element disagreement yields an `alt` copy
/// for the ledger's `discrepant` flag.
pub struct MergedRecord {
    pub best: CensusRecord,
    pub alt: Option<CensusRecord>,
}

pub fn merge_census(sbdb: Vec<CensusRecord>, mpc: Vec<CensusRecord>) -> Vec<MergedRecord> {
    let mut out: Vec<MergedRecord> = sbdb
        .into_iter()
        .map(|r| MergedRecord { best: r, alt: None })
        .collect();
    for m in mpc {
        let hit = out
            .iter_mut()
            .find(|o| o.best.aliases.iter().any(|a| m.aliases.contains(a)));
        match hit {
            Some(o) => {
                for a in &m.aliases {
                    if !o.best.aliases.contains(a) {
                        o.best.aliases.push(a.clone());
                    }
                }
                if o.best.n_oppositions.is_none() {
                    o.best.n_oppositions = m.n_oppositions;
                }
                if o.best.condition_code.is_none() {
                    o.best.condition_code = m.condition_code;
                }
                if o.best.arc_days.is_none() {
                    o.best.arc_days = m.arc_days;
                }
                if elements_discrepant(&o.best, &m) {
                    o.alt = Some(m);
                }
            }
            None => out.push(MergedRecord { best: m, alt: None }),
        }
    }
    out
}

/// Material disagreement: the brown2017-style drift allowances from the
/// refresh guard (a: 1% or 1 AU; q: 2%; angles 0.5°).
pub fn elements_discrepant(x: &CensusRecord, y: &CensusRecord) -> bool {
    let a_tol = (0.01 * x.a_au.abs()).max(1.0);
    (x.a_au - y.a_au).abs() > a_tol
        || (x.q_au - y.q_au).abs() > 0.02 * x.q_au.abs()
        || wrap_diff_deg(x.omega_deg, y.omega_deg) > 0.5
        || wrap_diff_deg(x.omega_big_deg, y.omega_big_deg) > 0.5
}

fn wrap_diff_deg(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d)
}

// ------------------------------------------------------------------ network

/// Fetch + gunzip + parse the MPC distant file. Returns (records,
/// last_modified header).
pub fn fetch_mpc_distant() -> Result<(Vec<CensusRecord>, Option<String>), String> {
    use std::io::Read;
    let mut resp = ureq::get(MPC_DISTANT_URL)
        .call()
        .map_err(|e| format!("mpc fetch: {e}"))?;
    let last_modified = resp
        .headers()
        .get("Last-Modified")
        .and_then(|v| v.to_str().ok())
        .map(str::to_string);
    let gz = resp
        .body_mut()
        .with_config()
        .limit(256 * 1024 * 1024)
        .read_to_vec()
        .map_err(|e| format!("mpc read: {e}"))?;
    let mut json = String::new();
    flate2::read::GzDecoder::new(&gz[..])
        .read_to_string(&mut json)
        .map_err(|e| format!("mpc gunzip: {e}"))?;
    Ok((parse_mpc_distant(&json)?, last_modified))
}

/// Fetch + parse the SBDB census.
pub fn fetch_sbdb_census() -> Result<Vec<CensusRecord>, String> {
    let mut resp = ureq::get(SBDB_QUERY_URL)
        .query("fields", SBDB_FIELDS)
        .query("sb-kind", SBDB_KIND)
        .query("sb-cdata", SBDB_CDATA)
        .call()
        .map_err(|e| format!("sbdb fetch: {e}"))?;
    let json = resp
        .body_mut()
        .with_config()
        .limit(64 * 1024 * 1024)
        .read_to_string()
        .map_err(|e| format!("sbdb read: {e}"))?;
    parse_sbdb_query(&json)
}
