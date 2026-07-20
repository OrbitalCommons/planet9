//! Layer-4 slow-mover linker (design/05): distance-hypothesis linking of
//! unassociated DIA detections into candidates moving like a body at
//! d ∈ [300, 1000] AU.
//!
//! The transform is heliolinc-style: assume a heliocentric distance r,
//! invert each topocentric detection to a heliocentric ecliptic direction
//! (λ, β); a real body at that r collapses to a slowly drifting point
//! (own motion ≤ 0.7″/day at 300 AU), while stars and wrong-r movers smear.
//! Clusters with ≥ 3 distinct nights are refit with the full 5-parameter
//! model (λ₀, β₀, 1/d, dλ/dt, dβ/dt) against the raw astrometry, must beat
//! a static-source model by Δχ² ≥ 25, and pass a loose photometric
//! plausibility gate before becoming candidate records (design/07).
//!
//! Approximations (documented in design/05): light-time retardation (a
//! near-constant few-day epoch shift at these distances) and differential
//! aberration (absorbed by Gaia-calibrated astrometry to first order) are
//! neglected; the injection harness uses the exact inverse of the same
//! transform, so measured completeness is internally consistent.

use p9_core::coords::observer::{EarthProvider, Time, Timescale};
use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, equatorial_to_ecliptic};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub const R_MIN_AU: f64 = 300.0;
pub const R_MAX_AU: f64 = 1000.0;
/// 1/d hypothesis step (AU⁻¹): keeps the wrong-hypothesis smear below ~2″
/// over a 90-day window (baseline ≤ 2 AU; design/05 grid math).
pub const INV_D_STEP: f64 = 4.0e-6;
/// Pair-linking tolerance (arcsec): transform residual + astrometric noise.
pub const LINK_TOL_ARCSEC: f64 = 3.0;
/// Maximum heliocentric drift rate allowed when pairing (arcsec/day): own
/// mean motion is ≤ 0.7″/day at 300 AU; 1.5 leaves margin without letting
/// asteroids in.
pub const DRIFT_MAX_ARCSEC_DAY: f64 = 1.5;
/// Astrometric σ per axis (arcsec) when the packet carries none.
pub const DEFAULT_SIGMA_ARCSEC: f64 = 0.15;
/// Candidate gates.
pub const MIN_NIGHTS: usize = 3;
pub const MAX_CHI2_PER_DOF: f64 = 3.0;
pub const MIN_STATIC_DELTA_CHI2: f64 = 25.0;
/// Loose implied-H plausibility window (mag): rejects fits whose distance
/// and brightness together demand an absurd body.
pub const H_MIN: f64 = -9.0;
pub const H_MAX: f64 = 5.0;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Detection {
    pub alert_id: u64,
    pub mjd: f64,
    pub ra: f64,
    pub dec: f64,
    pub band: String,
    pub psf_mag: Option<f64>,
    pub sigma_arcsec: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisitInfo {
    pub mjd: f64,
    pub band: String,
    pub ra: f64,
    pub dec: f64,
    pub radius_deg: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinkerInput {
    pub tile: u64,
    pub mjd_start: f64,
    pub mjd_end: f64,
    pub detections: Vec<Detection>,
    pub visits: Vec<VisitInfo>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CandidateFit {
    pub lambda0_deg: f64,
    pub beta0_deg: f64,
    pub inv_d_au: f64,
    pub d_au: f64,
    pub dlam_dt_deg_day: f64,
    pub dbeta_dt_deg_day: f64,
    pub chi2: f64,
    pub dof: usize,
    pub arc_days: f64,
    pub n_nights: usize,
    pub static_delta_chi2: f64,
    pub implied_h_mag: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Candidate {
    pub members: Vec<u64>,
    pub fit: CandidateFit,
}

/// Precomputed per-detection state shared by every hypothesis.
struct Prepared {
    /// Ecliptic-frame topocentric unit vector.
    u: [f64; 3],
    /// Earth heliocentric ecliptic position (AU).
    earth: [f64; 3],
    mjd: f64,
    sigma: f64,
    night: i64,
    idx: usize,
}

fn unit_from_ecliptic(lon: f64, lat: f64) -> [f64; 3] {
    let (cl, sl) = (lon.cos(), lon.sin());
    let (cb, sb) = (lat.cos(), lat.sin());
    [cb * cl, cb * sl, sb]
}

/// Chile-night grouping: the observing night spans the UTC-date boundary.
pub fn night_of(mjd: f64) -> i64 {
    (mjd - 0.5).floor() as i64
}

fn mjd_to_time(ts: &Timescale, mjd: f64) -> Time {
    // Alert epochs are TAI; TT − TAI = 32.184 s ≈ 3.7e-4 day — far below the
    // astrometric budget at these rates.
    ts.tt_jd(mjd + 2_400_000.5, None)
}

fn prepare(input: &LinkerInput, earth: &mut impl EarthProvider, ts: &Timescale) -> Vec<Prepared> {
    // One Earth state per unique epoch (visits share epochs across many
    // detections).
    let mut states: HashMap<i64, [f64; 3]> = HashMap::new();
    let mut out = Vec::with_capacity(input.detections.len());
    for (idx, det) in input.detections.iter().enumerate() {
        let key = (det.mjd * 86_400.0).round() as i64;
        let e = *states.entry(key).or_insert_with(|| {
            let p = earth.earth_position(&mjd_to_time(ts, det.mjd));
            [p.x, p.y, p.z]
        });
        let (lon, lat) = equatorial_to_ecliptic(det.ra.to_radians(), det.dec.to_radians());
        out.push(Prepared {
            u: unit_from_ecliptic(lon, lat),
            earth: e,
            mjd: det.mjd,
            sigma: det.sigma_arcsec.unwrap_or(DEFAULT_SIGMA_ARCSEC).max(0.02),
            night: night_of(det.mjd),
            idx,
        });
    }
    out
}

/// Heliocentric ecliptic (λ, β) of a detection under hypothesis distance
/// `r_au`. Returns None if the geometry is impossible (never for r > 1 AU).
fn helio_lonlat(p: &Prepared, r_au: f64) -> Option<(f64, f64)> {
    let e_dot_u = p.earth[0] * p.u[0] + p.earth[1] * p.u[1] + p.earth[2] * p.u[2];
    let e2 = p.earth[0] * p.earth[0] + p.earth[1] * p.earth[1] + p.earth[2] * p.earth[2];
    let disc = e_dot_u * e_dot_u + r_au * r_au - e2;
    if disc <= 0.0 {
        return None;
    }
    let rho = -e_dot_u + disc.sqrt();
    let x = [
        p.earth[0] + rho * p.u[0],
        p.earth[1] + rho * p.u[1],
        p.earth[2] + rho * p.u[2],
    ];
    let lon = x[1].atan2(x[0]);
    let lat = (x[2] / r_au).clamp(-1.0, 1.0).asin();
    Some((lon, lat))
}

/// Forward model: topocentric (ra, dec) in degrees of a body at heliocentric
/// (λ(t), β(t), r) seen from Earth position `e`.
fn forward_radec(lon: f64, lat: f64, r_au: f64, e: &[f64; 3]) -> (f64, f64) {
    let x = unit_from_ecliptic(lon, lat);
    let v = nalgebra::Vector3::new(r_au * x[0] - e[0], r_au * x[1] - e[1], r_au * x[2] - e[2]);
    ecliptic_vec_to_equatorial_deg(&v.normalize())
}

/// Union-find for clustering.
struct Dsu(Vec<usize>);
impl Dsu {
    fn new(n: usize) -> Self {
        Dsu((0..n).collect())
    }
    fn find(&mut self, i: usize) -> usize {
        if self.0[i] != i {
            let r = self.find(self.0[i]);
            self.0[i] = r;
        }
        self.0[i]
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0[ra] = rb;
        }
    }
}

/// One hypothesis pass: transform, drift-aware pair linking (λ-sorted sweep),
/// connected components with ≥ MIN_NIGHTS distinct nights.
fn cluster_at_hypothesis(prep: &[Prepared], r_au: f64) -> Vec<Vec<usize>> {
    const ARCSEC: f64 = 1.0 / 206_265.0;
    let tol = LINK_TOL_ARCSEC * ARCSEC;
    let drift = DRIFT_MAX_ARCSEC_DAY * ARCSEC;

    let mut pts: Vec<(f64, f64, f64, usize)> = Vec::with_capacity(prep.len());
    for p in prep {
        if let Some((lon, lat)) = helio_lonlat(p, r_au) {
            pts.push((lon, lat, p.mjd, p.idx));
        }
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let span = pts.last().map(|l| l.2).unwrap_or(0.0) - pts.first().map(|f| f.2).unwrap_or(0.0);
    let max_sep = tol + drift * span.max(1.0);

    let mut dsu = Dsu::new(prep.len());
    for i in 0..pts.len() {
        let (li, bi, ti, idx_i) = pts[i];
        let cos_b = bi.cos().max(1e-6);
        for &(lj, bj, tj, idx_j) in pts.iter().skip(i + 1) {
            if (lj - li) * cos_b > max_sep {
                break; // λ-sorted sweep bound
            }
            let dt = (tj - ti).abs();
            let allowed = tol + drift * dt.max(0.2);
            let dlam = (lj - li) * cos_b;
            let dbet = bj - bi;
            if dlam * dlam + dbet * dbet <= allowed * allowed {
                dsu.union(idx_i, idx_j);
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for p in prep {
        groups.entry(dsu.find(p.idx)).or_default().push(p.idx);
    }
    groups
        .into_values()
        .filter(|g| {
            let nights: std::collections::HashSet<i64> = g.iter().map(|&i| prep[i].night).collect();
            nights.len() >= MIN_NIGHTS
        })
        .collect()
}

/// 5-parameter fit for one member set: full-range 1/d search (coarse scan
/// then fine refinement — the clustering stage only proposes member sets;
/// residual-parallax chains can link a true object at badly wrong
/// hypotheses, so the fit must not trust the seed), linear (λ, β) drift
/// regression in the transformed frame, χ² evaluated by forward-modeling
/// back to the raw astrometry.
fn fit_members(prep: &[Prepared], members: &[usize]) -> Option<(CandidateFit, f64)> {
    let t0 = members.iter().map(|&i| prep[i].mjd).sum::<f64>() / members.len() as f64;

    let eval = |inv_d: f64| -> Option<(f64, [f64; 4])> {
        let r = 1.0 / inv_d;
        // Transform members; regress λ, β on Δt.
        let mut lam = Vec::new();
        let mut bet = Vec::new();
        let mut dts = Vec::new();
        let mut lam_ref = None;
        for &i in members {
            let (l, b) = helio_lonlat(&prep[i], r)?;
            let lr = *lam_ref.get_or_insert(l);
            let mut l = l;
            while l - lr > std::f64::consts::PI {
                l -= std::f64::consts::TAU;
            }
            while l - lr < -std::f64::consts::PI {
                l += std::f64::consts::TAU;
            }
            lam.push(l);
            bet.push(b);
            dts.push(prep[i].mjd - t0);
        }
        let reg = |ys: &[f64]| -> (f64, f64) {
            let n = ys.len() as f64;
            let st: f64 = dts.iter().sum();
            let stt: f64 = dts.iter().map(|t| t * t).sum();
            let sy: f64 = ys.iter().sum();
            let sty: f64 = dts.iter().zip(ys).map(|(t, y)| t * y).sum();
            let denom = n * stt - st * st;
            if denom.abs() < 1e-12 {
                (sy / n, 0.0)
            } else {
                ((sy * stt - st * sty) / denom, (n * sty - st * sy) / denom)
            }
        };
        let (l0, dl) = reg(&lam);
        let (b0, db) = reg(&bet);
        // χ² against the raw astrometry.
        let mut chi2 = 0.0;
        for (k, &i) in members.iter().enumerate() {
            let p = &prep[i];
            let (ra_p, dec_p) = forward_radec(l0 + dl * dts[k], b0 + db * dts[k], r, &p.earth);
            let det = /* observed */ (p.u, p.mjd);
            let _ = det;
            let (ra_o, dec_o) = {
                // Recover observed ra/dec from the stored unit vector.
                let v = nalgebra::Vector3::new(p.u[0], p.u[1], p.u[2]);
                ecliptic_vec_to_equatorial_deg(&v)
            };
            let cosd = dec_o.to_radians().cos();
            let mut dra = (ra_p - ra_o).rem_euclid(360.0);
            if dra > 180.0 {
                dra -= 360.0;
            }
            let dra_as = dra * cosd * 3600.0;
            let ddec_as = (dec_p - dec_o) * 3600.0;
            let s = p.sigma;
            chi2 += (dra_as / s).powi(2) + (ddec_as / s).powi(2);
        }
        Some((chi2, [l0, b0, dl, db]))
    };

    // Stage 1: coarse scan of the full 1/d range; stage 2: fine refinement
    // (INV_D_STEP resolution) around the coarse minimum.
    let (inv_lo, inv_hi) = (1.0 / R_MAX_AU, 1.0 / R_MIN_AU);
    let coarse_step = (inv_hi - inv_lo) / 60.0;
    let mut best: Option<(f64, f64, [f64; 4])> = None;
    let mut consider = |inv_d: f64, best: &mut Option<(f64, f64, [f64; 4])>| {
        if let Some((chi2, params)) = eval(inv_d) {
            if best.as_ref().is_none_or(|(bc, _, _)| chi2 < *bc) {
                *best = Some((chi2, inv_d, params));
            }
        }
    };
    for k in 0..=60 {
        consider(inv_lo + k as f64 * coarse_step, &mut best);
    }
    let coarse_center = best.as_ref()?.1;
    let fine_span = (coarse_step / INV_D_STEP).ceil() as i64;
    for k in -fine_span..=fine_span {
        let inv_d = (coarse_center + k as f64 * INV_D_STEP).clamp(inv_lo, inv_hi);
        consider(inv_d, &mut best);
    }
    let (chi2, inv_d, [l0, b0, dl, db]) = best?;

    // Static-source model: best fixed topocentric position.
    let (mut sra, mut sdec, mut n) = (0.0, 0.0, 0.0);
    let obs: Vec<(f64, f64)> = members
        .iter()
        .map(|&i| {
            let v = nalgebra::Vector3::new(prep[i].u[0], prep[i].u[1], prep[i].u[2]);
            ecliptic_vec_to_equatorial_deg(&v)
        })
        .collect();
    // Handle RA wrap by vector mean.
    let (mut x, mut y) = (0.0, 0.0);
    for &(ra, dec) in &obs {
        x += ra.to_radians().cos();
        y += ra.to_radians().sin();
        sdec += dec;
        n += 1.0;
    }
    sra = y.atan2(x).to_degrees().rem_euclid(360.0);
    sdec /= n;
    let mut chi2_static = 0.0;
    for (k, &i) in members.iter().enumerate() {
        let (ra_o, dec_o) = obs[k];
        let cosd = dec_o.to_radians().cos();
        let mut dra = (sra - ra_o).rem_euclid(360.0);
        if dra > 180.0 {
            dra -= 360.0;
        }
        let s = prep[i].sigma;
        chi2_static += (dra * cosd * 3600.0 / s).powi(2) + ((sdec - dec_o) * 3600.0 / s).powi(2);
    }

    let mjds: Vec<f64> = members.iter().map(|&i| prep[i].mjd).collect();
    let arc = mjds.iter().cloned().fold(f64::MIN, f64::max)
        - mjds.iter().cloned().fold(f64::MAX, f64::min);
    let nights: std::collections::HashSet<i64> = members.iter().map(|&i| prep[i].night).collect();
    let dof = 2 * members.len() - 5;

    Some((
        CandidateFit {
            lambda0_deg: l0.to_degrees().rem_euclid(360.0),
            beta0_deg: b0.to_degrees(),
            inv_d_au: inv_d,
            d_au: 1.0 / inv_d,
            dlam_dt_deg_day: dl.to_degrees(),
            dbeta_dt_deg_day: db.to_degrees(),
            chi2,
            dof,
            arc_days: arc,
            n_nights: nights.len(),
            static_delta_chi2: chi2_static - chi2,
            implied_h_mag: None,
        },
        chi2_static,
    ))
}

/// Run the linker over one tile-window input. `detections` should already be
/// the unassociated selection; the caller decides that policy.
pub fn link(input: &LinkerInput, earth: &mut impl EarthProvider, ts: &Timescale) -> Vec<Candidate> {
    let prep = prepare(input, earth, ts);
    if prep.len() < MIN_NIGHTS {
        return Vec::new();
    }

    let inv_lo = 1.0 / R_MAX_AU;
    let inv_hi = 1.0 / R_MIN_AU;
    let steps = ((inv_hi - inv_lo) / INV_D_STEP).ceil() as usize;

    // Collect raw (member-set → seed) candidates across hypotheses, deduped
    // by member overlap.
    let mut raw: Vec<(Vec<usize>, f64)> = Vec::new();
    for k in 0..=steps {
        let inv_d = inv_lo + k as f64 * INV_D_STEP;
        let groups = cluster_at_hypothesis(&prep, 1.0 / inv_d);
        if !groups.is_empty() && std::env::var("P9_LINK_DEBUG").is_ok() {
            eprintln!(
                "k={k} inv_d={inv_d:.6e} r={:.0}: {} groups",
                1.0 / inv_d,
                groups.len()
            );
        }
        for mut group in groups {
            group.sort_unstable();
            match raw.iter_mut().find(|(g, _)| overlap_frac(g, &group) >= 0.5) {
                Some((g, _)) => {
                    if group.len() > g.len() {
                        *g = group;
                    }
                }
                None => raw.push((group, inv_d)),
            }
        }
    }

    let mut out = Vec::new();
    for (members, _seed) in raw {
        let Some((mut fit, _)) = fit_members(&prep, &members) else {
            continue;
        };
        if fit.dof == 0 || fit.chi2 / fit.dof as f64 > MAX_CHI2_PER_DOF {
            continue;
        }
        if fit.static_delta_chi2 < MIN_STATIC_DELTA_CHI2 {
            continue;
        }
        // Photometric plausibility: implied absolute magnitude at the fitted
        // distance from the mean detected magnitude (opposition law shape;
        // loose window).
        let mags: Vec<f64> = members
            .iter()
            .filter_map(|&i| input.detections[i].psf_mag)
            .collect();
        if !mags.is_empty() {
            let mean_m = mags.iter().sum::<f64>() / mags.len() as f64;
            let d = fit.d_au;
            let h = mean_m - 5.0 * (d * (d - 1.0)).log10();
            fit.implied_h_mag = Some(h);
            if !(H_MIN..=H_MAX).contains(&h) {
                continue;
            }
        }
        out.push(Candidate {
            members: members
                .iter()
                .map(|&i| input.detections[i].alert_id)
                .collect(),
            fit,
        });
    }
    out.sort_by(|a, b| a.fit.chi2.partial_cmp(&b.fit.chi2).unwrap());
    out
}

fn overlap_frac(a: &[usize], b: &[usize]) -> f64 {
    let sa: std::collections::HashSet<_> = a.iter().collect();
    let inter = b.iter().filter(|x| sa.contains(x)).count();
    inter as f64 / a.len().min(b.len()).max(1) as f64
}

// ------------------------------------------------------------- calibration

#[derive(Debug, Clone, Serialize)]
pub struct InjectionResult {
    pub injected: usize,
    pub recovered: usize,
    pub efficiency: f64,
    /// Per-injection detail: (d_au, v_mag, n_detections, recovered,
    /// fitted_d_if_recovered).
    pub detail: Vec<(f64, f64, usize, bool, Option<f64>)>,
}

/// Injection/recovery harness (SR-6): synthetic bodies at d ∈ [300, 1000] AU
/// dropped onto the input's REAL visit list, detected through the per-band
/// logistic efficiency, then recovered (or not) by the full linker. The
/// forward model is the exact inverse of the linking transform, so this
/// measures pipeline completeness, not transform truth.
pub fn calibrate(
    input: &LinkerInput,
    earth: &mut impl EarthProvider,
    ts: &Timescale,
    n_inject: usize,
    seed: u64,
) -> InjectionResult {
    use p9_core::analysis::surveys::logistic_efficiency;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let fiducial = |band: &str| match band {
        "u" => 23.9,
        "g" => 24.8,
        "r" => 24.3,
        "i" => 23.9,
        "z" => 23.3,
        _ => 22.1,
    };

    let mut detail = Vec::new();
    let mut recovered_count = 0usize;
    for inj in 0..n_inject {
        // Truth: distance, position near a random visit center, own motion
        // bounded by the circular mean motion at d.
        let d = rng.gen_range(R_MIN_AU..R_MAX_AU);
        let v_mag = rng.gen_range(20.0..23.5);
        let vref = &input.visits[rng.gen_range(0..input.visits.len())];
        let t_ref = mjd_to_time(ts, vref.mjd);
        let e_ref = earth.earth_position(&t_ref);
        // Place truth AT the reference visit center (transformed).
        let (lon_t, lat_t) = equatorial_to_ecliptic(vref.ra.to_radians(), vref.dec.to_radians());
        let p_ref = Prepared {
            u: unit_from_ecliptic(lon_t, lat_t),
            earth: [e_ref.x, e_ref.y, e_ref.z],
            mjd: vref.mjd,
            sigma: DEFAULT_SIGMA_ARCSEC,
            night: night_of(vref.mjd),
            idx: 0,
        };
        let Some((l0, b0)) = helio_lonlat(&p_ref, d) else {
            continue;
        };
        let n_mean = 0.9856 / d.powf(1.5); // deg/day circular mean motion
        let dl = rng.gen_range(-n_mean..n_mean);
        let db = rng.gen_range(-n_mean / 2.0..n_mean / 2.0);

        // Generate detections at every visit whose cone contains the body.
        let mut synth = Vec::new();
        for v in &input.visits {
            let t = mjd_to_time(ts, v.mjd);
            let e = earth.earth_position(&t);
            let dt = v.mjd - vref.mjd;
            let (ra_p, dec_p) = forward_radec(
                l0 + (dl * dt).to_radians(),
                b0 + (db * dt).to_radians(),
                d,
                &[e.x, e.y, e.z],
            );
            let sep = angular_sep_deg(ra_p, dec_p, v.ra, v.dec);
            if sep > v.radius_deg {
                continue;
            }
            let eps = logistic_efficiency(v_mag, fiducial(&v.band), 4.0);
            if rng.gen::<f64>() >= eps {
                continue;
            }
            let noise = DEFAULT_SIGMA_ARCSEC / 3600.0;
            synth.push(Detection {
                alert_id: 1_000_000_000 + (inj * 1000 + synth.len()) as u64,
                mjd: v.mjd,
                ra: (ra_p + rng.gen_range(-noise..noise) / dec_p.to_radians().cos().max(1e-6))
                    .rem_euclid(360.0),
                dec: dec_p + rng.gen_range(-noise..noise),
                band: v.band.clone(),
                psf_mag: Some(v_mag),
                sigma_arcsec: Some(DEFAULT_SIGMA_ARCSEC),
            });
        }
        let nights: std::collections::HashSet<i64> =
            synth.iter().map(|s| night_of(s.mjd)).collect();
        let injectable = nights.len() >= MIN_NIGHTS;

        // Run the linker on background + synthetic.
        let mut merged = input.clone();
        let synth_ids: std::collections::HashSet<u64> = synth.iter().map(|s| s.alert_id).collect();
        let n_synth = synth.len();
        merged.detections.extend(synth);
        let cands = link(&merged, earth, ts);
        let hit = cands.iter().find(|c| {
            let got = c.members.iter().filter(|m| synth_ids.contains(m)).count();
            got * 2 >= n_synth.max(1) && (c.fit.d_au / d - 1.0).abs() < 0.15
        });
        let rec = injectable && hit.is_some();
        if rec {
            recovered_count += 1;
        }
        detail.push((d, v_mag, n_synth, rec, hit.map(|c| c.fit.d_au)));
    }
    let injectable_n = detail
        .iter()
        .filter(|(_, _, n, _, _)| {
            *n > 0 // at least detected somewhere
        })
        .count();
    let _ = injectable_n;
    let injected = detail.len();
    InjectionResult {
        injected,
        recovered: recovered_count,
        efficiency: recovered_count as f64 / injected.max(1) as f64,
        detail,
    }
}

fn angular_sep_deg(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let (r1, d1, r2, d2) = (
        ra1.to_radians(),
        dec1.to_radians(),
        ra2.to_radians(),
        dec2.to_radians(),
    );
    let c = d1.sin() * d2.sin() + d1.cos() * d2.cos() * (r1 - r2).cos();
    c.clamp(-1.0, 1.0).acos().to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::coords::observer::CircularEarth;

    /// A dense synthetic visit list: same field revisited over 5 nights in a
    /// 60-day window (r band), generous cone.
    fn test_visits() -> Vec<VisitInfo> {
        let mut v = Vec::new();
        for k in 0..10 {
            v.push(VisitInfo {
                mjd: 61200.25 + (k as f64) * 6.0,
                band: "r".into(),
                ra: 40.0,
                dec: -10.0,
                radius_deg: 1.75,
            });
        }
        v
    }

    fn empty_input() -> LinkerInput {
        LinkerInput {
            tile: 0,
            mjd_start: 61200.0,
            mjd_end: 61260.0,
            detections: Vec::new(),
            visits: test_visits(),
        }
    }

    /// Forward-generate a body at known distance and check the full linker
    /// recovers it with the right distance.
    #[test]
    fn recovers_synthetic_slow_mover() {
        let mut earth = CircularEarth;
        let ts = Timescale::default();
        let d_true = 600.0;

        let mut input = empty_input();
        let visits = input.visits.clone();
        // Truth placed at the field center via the same transform.
        let t0 = mjd_to_time(&ts, visits[0].mjd);
        let e0 = earth.earth_position(&t0);
        let (lon, lat) =
            equatorial_to_ecliptic(visits[0].ra.to_radians(), visits[0].dec.to_radians());
        let p0 = Prepared {
            u: unit_from_ecliptic(lon, lat),
            earth: [e0.x, e0.y, e0.z],
            mjd: visits[0].mjd,
            sigma: DEFAULT_SIGMA_ARCSEC,
            night: night_of(visits[0].mjd),
            idx: 0,
        };
        let (l0, b0) = helio_lonlat(&p0, d_true).unwrap();
        let n_mean = (0.9856 / d_true.powf(1.5)).to_radians(); // rad/day
        for (k, v) in visits.iter().enumerate() {
            let t = mjd_to_time(&ts, v.mjd);
            let e = earth.earth_position(&t);
            let dt = v.mjd - visits[0].mjd;
            let (ra, dec) = forward_radec(l0 + 0.5 * n_mean * dt, b0, d_true, &[e.x, e.y, e.z]);
            input.detections.push(Detection {
                alert_id: k as u64,
                mjd: v.mjd,
                ra,
                dec,
                band: "r".into(),
                psf_mag: Some(21.8),
                sigma_arcsec: Some(0.1),
            });
        }

        let cands = link(&input, &mut earth, &ts);
        assert_eq!(cands.len(), 1, "exactly one candidate");
        let c = &cands[0];
        assert_eq!(c.members.len(), 10);
        assert!(
            (c.fit.d_au / d_true - 1.0).abs() < 0.10,
            "fitted d = {:.0}",
            c.fit.d_au
        );
        assert!(c.fit.static_delta_chi2 >= MIN_STATIC_DELTA_CHI2);
        assert!(c.fit.chi2 / c.fit.dof as f64 <= MAX_CHI2_PER_DOF);
        // Implied H at 600 AU and V 21.8 is planet-like.
        let h = c.fit.implied_h_mag.unwrap();
        assert!((-8.0..0.0).contains(&h), "H = {h:.1}");
    }

    #[test]
    fn static_source_is_rejected() {
        let mut earth = CircularEarth;
        let ts = Timescale::default();
        let mut input = empty_input();
        for (k, v) in input.visits.clone().iter().enumerate() {
            input.detections.push(Detection {
                alert_id: k as u64,
                mjd: v.mjd,
                ra: 40.123,
                dec: -10.234, // perfectly static: a star/AGN
                band: "r".into(),
                psf_mag: Some(21.0),
                sigma_arcsec: Some(0.1),
            });
        }
        let cands = link(&input, &mut earth, &ts);
        assert!(
            cands.is_empty(),
            "static source must not survive: {:?}",
            cands[0].fit
        );
    }

    #[test]
    fn two_nights_are_not_enough() {
        let mut earth = CircularEarth;
        let ts = Timescale::default();
        let mut input = empty_input();
        input.visits.truncate(2);
        let visits = input.visits.clone();
        // Perfect mover but only 2 nights.
        let t0 = mjd_to_time(&ts, visits[0].mjd);
        let e0 = earth.earth_position(&t0);
        let (lon, lat) =
            equatorial_to_ecliptic(visits[0].ra.to_radians(), visits[0].dec.to_radians());
        let p0 = Prepared {
            u: unit_from_ecliptic(lon, lat),
            earth: [e0.x, e0.y, e0.z],
            mjd: visits[0].mjd,
            sigma: 0.1,
            night: night_of(visits[0].mjd),
            idx: 0,
        };
        let (l0, b0) = helio_lonlat(&p0, 500.0).unwrap();
        for (k, v) in visits.iter().enumerate() {
            let t = mjd_to_time(&ts, v.mjd);
            let e = earth.earth_position(&t);
            let (ra, dec) = forward_radec(l0, b0, 500.0, &[e.x, e.y, e.z]);
            input.detections.push(Detection {
                alert_id: k as u64,
                mjd: v.mjd,
                ra,
                dec,
                band: "r".into(),
                psf_mag: Some(22.0),
                sigma_arcsec: Some(0.1),
            });
        }
        assert!(link(&input, &mut earth, &ts).is_empty());
    }

    #[test]
    fn nearby_decoy_outside_distance_range_is_rejected() {
        // A 120 AU body moves too fast in the transformed frame at every
        // hypothesis in [300, 1000]: residual parallax smears it beyond the
        // drift-aware linking radius.
        let mut earth = CircularEarth;
        let ts = Timescale::default();
        let d_decoy = 120.0;
        let mut input = empty_input();
        let visits = input.visits.clone();
        let t0 = mjd_to_time(&ts, visits[0].mjd);
        let e0 = earth.earth_position(&t0);
        let (lon, lat) =
            equatorial_to_ecliptic(visits[0].ra.to_radians(), visits[0].dec.to_radians());
        let p0 = Prepared {
            u: unit_from_ecliptic(lon, lat),
            earth: [e0.x, e0.y, e0.z],
            mjd: visits[0].mjd,
            sigma: 0.1,
            night: night_of(visits[0].mjd),
            idx: 0,
        };
        let (l0, b0) = helio_lonlat(&p0, d_decoy).unwrap();
        for (k, v) in visits.iter().enumerate() {
            let t = mjd_to_time(&ts, v.mjd);
            let e = earth.earth_position(&t);
            let (ra, dec) = forward_radec(l0, b0, d_decoy, &[e.x, e.y, e.z]);
            input.detections.push(Detection {
                alert_id: k as u64,
                mjd: v.mjd,
                ra,
                dec,
                band: "r".into(),
                psf_mag: Some(21.0),
                sigma_arcsec: Some(0.1),
            });
        }
        let cands = link(&input, &mut earth, &ts);
        assert!(
            cands.is_empty(),
            "120 AU decoy must not fit in [300, 1000]: {:?}",
            cands[0].fit
        );
    }

    #[test]
    fn calibration_meets_sr6_on_dense_cadence() {
        let mut earth = CircularEarth;
        let ts = Timescale::default();
        let input = empty_input();
        let res = calibrate(&input, &mut earth, &ts, 40, 7);
        // SR-6: >= 90% for V <= 23 with >= 3 covered nights. The dense test
        // cadence guarantees coverage, so overall efficiency should clear
        // the gate for the V <= 23 subset.
        let sub: Vec<_> = res
            .detail
            .iter()
            .filter(|(_, v, n, _, _)| *v <= 23.0 && *n >= MIN_NIGHTS)
            .collect();
        let rec = sub.iter().filter(|(_, _, _, r, _)| *r).count();
        assert!(
            !sub.is_empty(),
            "test cadence should yield linkable injections"
        );
        let eff = rec as f64 / sub.len() as f64;
        assert!(
            eff >= 0.9,
            "SR-6: efficiency {eff:.2} on {} linkable",
            sub.len()
        );
    }
}
