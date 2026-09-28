//! Film export for `p9-2016-resonance-prediction`: the numbers its scene and ledger entry draw.

use p9_2016_resonance_prediction::{
    DEFAULT_P_MAX, DEFAULT_Q_MAX, MALHOTRA_2016_A9_AU, MALHOTRA_2016_PERIOD_YR,
    MILLHOLLAND_2016_A9_AU, best_fit_a9, malhotra_constraining_set, nearest_resonance,
    orbital_period_years, scan_a9,
};
use serde_json::{Value, json};

/// Scan range and resolution for the candidate Planet Nine semi-major axis.
const A_MIN: f64 = 400.0;
const A_MAX: f64 = 900.0;
const FIT_STEPS: usize = 5001;
const CURVE_STEPS: usize = 501;

fn matches(a9: f64) -> Vec<Value> {
    malhotra_constraining_set()
        .iter()
        .map(|e| {
            let m = nearest_resonance(e.a, a9, DEFAULT_P_MAX, DEFAULT_Q_MAX)
                .expect("every constraining ETNO lies inside the scanned a9");
            json!({
                "name": e.name,
                "a_au": e.a,
                "period_yr": orbital_period_years(e.a),
                "p": m.p,
                "q": m.q,
                "period_ratio": m.period_ratio,
                "residual": m.residual,
                "a_res_au": a9 * (m.q as f64 / m.p as f64).powf(2.0 / 3.0),
            })
        })
        .collect()
}

pub fn export() -> Value {
    let set: Vec<f64> = malhotra_constraining_set().iter().map(|e| e.a).collect();
    let best = best_fit_a9(&set, A_MIN, A_MAX, FIT_STEPS, DEFAULT_P_MAX, DEFAULT_Q_MAX);

    let curve = scan_a9(
        &set,
        A_MIN,
        A_MAX,
        CURVE_STEPS,
        DEFAULT_P_MAX,
        DEFAULT_Q_MAX,
    );
    let finite: Vec<_> = curve.iter().filter(|s| s.residual.is_finite()).collect();

    // How much of the scanned range fits at least as well as the published
    // solution does.
    let fine = scan_a9(&set, A_MIN, A_MAX, FIT_STEPS, DEFAULT_P_MAX, DEFAULT_Q_MAX);
    let published_residual = fine
        .iter()
        .min_by(|x, y| {
            (x.a9_au - MALHOTRA_2016_A9_AU)
                .abs()
                .partial_cmp(&(y.a9_au - MALHOTRA_2016_A9_AU).abs())
                .unwrap()
        })
        .map(|s| s.residual)
        .unwrap();
    let n_finite = fine.iter().filter(|s| s.residual.is_finite()).count();
    let n_as_good = fine
        .iter()
        .filter(|s| s.residual <= published_residual)
        .count();

    // The N:1 and N:2 ladder of resonant semi-major axes for the best fit.
    let mut ladder = Vec::new();
    for q in 1..=DEFAULT_Q_MAX {
        for p in (q + 1)..=DEFAULT_P_MAX {
            if p % q == 0 && q > 1 {
                continue;
            }
            let a_res = best.a9_au * (q as f64 / p as f64).powf(2.0 / 3.0);
            ladder.push(json!({"p": p, "q": q, "a_res_au": a_res}));
        }
    }

    json!({
        "best_a9_au": best.a9_au,
        "best_period_yr": best.period_yr,
        "best_residual": best.residual,
        "published_a9_au": MALHOTRA_2016_A9_AU,
        "published_period_yr": MALHOTRA_2016_PERIOD_YR,
        "published_residual": published_residual,
        "millholland_a9_au": MILLHOLLAND_2016_A9_AU,
        "fraction_as_good_as_published": n_as_good as f64 / n_finite as f64,
        "scan": {
            "a9_au": finite.iter().map(|s| s.a9_au).collect::<Vec<_>>(),
            "residual": finite.iter().map(|s| s.residual).collect::<Vec<_>>(),
        },
        "matches_best": matches(best.a9_au),
        "matches_published": matches(MALHOTRA_2016_A9_AU),
        "ladder": ladder,
    })
}
