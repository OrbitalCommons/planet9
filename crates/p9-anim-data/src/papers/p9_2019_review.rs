//! Film export for `p9-2019-review`: the numbers its scene and ledger entry draw.

use p9_2019_clustering::kbo_sample::paper_sample_a230;
use p9_2019_review::detection_prospects::{brightness_table, survey_limits};
use p9_2019_review::parameter_survey::{ParameterGrid, critical_semi_major_axis};
use p9_2019_review::revised_parameters::{
    P9ParameterSet, ParameterRange, best_fit_5me, best_fit_10me, original_2016, revised_2019,
};
use p9_core::types::P9Params;
use serde_json::{Value, json};

use super::p9_2016_evidence::{orbit_json, p9_orbit_json};

fn range_json(r: &ParameterRange) -> Value {
    json!({"min": r.min, "max": r.max, "best": r.best})
}

fn set_json(s: &P9ParameterSet) -> Value {
    json!({
        "label": s.label,
        "mass": range_json(&s.mass_earth),
        "a": range_json(&s.a),
        "e": range_json(&s.e),
        "i": range_json(&s.i_deg),
    })
}

pub fn export() -> Value {
    let before = P9Params::nominal_2016();
    let after = P9Params::revised_2019();
    let a_c_before = critical_semi_major_axis(&before);
    let a_c_after = critical_semi_major_axis(&after);

    let fits: Vec<Value> = best_fit_5me()
        .iter()
        .chain(best_fit_10me().iter())
        .map(|p| {
            let mut orbit = p9_orbit_json("fit", p);
            orbit["a_c"] = json!(critical_semi_major_axis(p));
            orbit
        })
        .collect();

    let brightness: Vec<Value> = brightness_table()
        .iter()
        .map(|b| {
            json!({
                "label": b.label,
                "mass": b.mass_earth,
                "v_perihelion_bright": b.v_perihelion_bright,
                "v_perihelion_faint": b.v_perihelion_faint,
                "v_aphelion_bright": b.v_aphelion_bright,
                "v_aphelion_faint": b.v_aphelion_faint,
            })
        })
        .collect();
    let surveys: Vec<Value> = survey_limits()
        .iter()
        .map(|s| json!({"name": s.name, "v_limit": s.v_limit}))
        .collect();

    let objects: Vec<Value> = paper_sample_a230()
        .iter()
        .map(|k| orbit_json(k.name, &k.elements))
        .collect();

    json!({
        "mass": after.mass_earth,
        "a": after.a,
        "e": after.e,
        "i": after.i.to_degrees(),
        "orbit_2016": p9_orbit_json("Batygin & Brown 2016", &before),
        "orbit_2019": p9_orbit_json("Batygin et al. 2019", &after),
        "ranges_2016": set_json(&original_2016()),
        "ranges_2019": set_json(&revised_2019()),
        "a_c_2016": a_c_before,
        "a_c": a_c_after,
        "fits": fits,
        "n_simulations": ParameterGrid::paper_grid().count_viable(),
        "brightness": brightness,
        "surveys": surveys,
        "objects": objects,
    })
}
