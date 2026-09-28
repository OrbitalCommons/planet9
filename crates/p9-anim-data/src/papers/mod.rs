//! Per-paper film exports: one module per reproduction crate, each returning the
//! values its scene and ledger entry read from `anim.json -> papers -> <crate>`.

use std::collections::BTreeMap;

use serde_json::Value;

pub mod p9_2016_cassini_ranging;
pub mod p9_2016_commensurabilities;
pub mod p9_2016_constraints;
pub mod p9_2016_cowan_thermal;
pub mod p9_2016_evidence;
pub mod p9_2016_fortney_thermal;
pub mod p9_2016_holman_payne;
pub mod p9_2016_inclination_instability;
pub mod p9_2016_inclined_tnos;
pub mod p9_2016_iorio_precession;
pub mod p9_2016_linder_evolution;
pub mod p9_2016_obliquity;
pub mod p9_2016_obliquity_gomes;
pub mod p9_2016_obliquity_lai;
pub mod p9_2016_resonance_prediction;
pub mod p9_2016_secular_resonance;
pub mod p9_2016_sheppard_etnos;
pub mod p9_2016_wise_coadd;
pub mod p9_2017_bias;
pub mod p9_2017_dynamics;
pub mod p9_2017_ossos_bias;
pub mod p9_2017_resonance_hopping;
pub mod p9_2018_bp519;
pub mod p9_2018_chaotic_dynamics;
pub mod p9_2018_kuiper_belt;
pub mod p9_2018_low_perihelion;
pub mod p9_2018_resonance;
pub mod p9_2018_secular_dynamics;
pub mod p9_2018_wise_search;
pub mod p9_2019_clustering;
pub mod p9_2019_ossos_scattering;
pub mod p9_2019_review;
pub mod p9_2019_selfgrav_disk;
pub mod p9_2019_tess_yield;
pub mod p9_2020_clement_precession;
pub mod p9_2020_des_isotropy;
pub mod p9_2020_resonance_hopping;
pub mod p9_2020_secular_octupole;
pub mod p9_2020_tess_shiftstack;
pub mod p9_2021_act_mm;
pub mod p9_2021_des_catalog;
pub mod p9_2021_detached_inclinations;
pub mod p9_2021_napier_critique;
pub mod p9_2021_oort_cloud;
pub mod p9_2021_orbit;
pub mod p9_2021_perihelion_gap;
pub mod p9_2021_stability;
pub mod p9_2021_ztf;
pub mod p9_2022_des;
pub mod p9_2022_iras_candidate;
pub mod p9_2022_uranus_tilt;
pub mod p9_2023_lsst_strategy;
pub mod p9_2023_mond_efe;
pub mod p9_2024_neptune_crossing;
pub mod p9_2024_oort_selfgrav;
pub mod p9_2024_panstarrs;
pub mod p9_2024_primordial_alignment;
pub mod p9_2024_siraj_orbit;
pub mod p9_2025_akari_refutation;
pub mod p9_2025_clustering;
pub mod p9_2025_iras_akari;
pub mod p9_2025_new_discoveries;
pub mod p9_2025_parallax_search;
pub mod p9_2025_perturbation;
pub mod p9_2025_planet_y;
pub mod p9_2025_ps1_holman;
pub mod p9_2025_russell_albedo;
pub mod p9_2025_simons_forecast;
pub mod p9_2025_stacking;
pub mod p9_2025_stellar_flybys;
pub mod p9_2026_alpha_slope;
pub mod p9_2026_apsidal_clustering;
pub mod p9_2026_cluster_inclinations;
pub mod p9_2026_flyby_evolution;
pub mod p9_2026_iorio_precession;
pub mod p9_2026_stellar_companions;

/// Every paper's export, keyed by crate name.
pub fn all() -> BTreeMap<&'static str, Value> {
    BTreeMap::from([
        ("p9-2016-cassini-ranging", p9_2016_cassini_ranging::export()),
        (
            "p9-2016-commensurabilities",
            p9_2016_commensurabilities::export(),
        ),
        ("p9-2016-constraints", p9_2016_constraints::export()),
        ("p9-2016-cowan-thermal", p9_2016_cowan_thermal::export()),
        ("p9-2016-evidence", p9_2016_evidence::export()),
        ("p9-2016-fortney-thermal", p9_2016_fortney_thermal::export()),
        ("p9-2016-holman-payne", p9_2016_holman_payne::export()),
        (
            "p9-2016-inclination-instability",
            p9_2016_inclination_instability::export(),
        ),
        ("p9-2016-inclined-tnos", p9_2016_inclined_tnos::export()),
        (
            "p9-2016-iorio-precession",
            p9_2016_iorio_precession::export(),
        ),
        (
            "p9-2016-linder-evolution",
            p9_2016_linder_evolution::export(),
        ),
        ("p9-2016-obliquity", p9_2016_obliquity::export()),
        ("p9-2016-obliquity-gomes", p9_2016_obliquity_gomes::export()),
        ("p9-2016-obliquity-lai", p9_2016_obliquity_lai::export()),
        (
            "p9-2016-resonance-prediction",
            p9_2016_resonance_prediction::export(),
        ),
        (
            "p9-2016-secular-resonance",
            p9_2016_secular_resonance::export(),
        ),
        ("p9-2016-sheppard-etnos", p9_2016_sheppard_etnos::export()),
        ("p9-2016-wise-coadd", p9_2016_wise_coadd::export()),
        ("p9-2017-bias", p9_2017_bias::export()),
        ("p9-2017-dynamics", p9_2017_dynamics::export()),
        ("p9-2017-ossos-bias", p9_2017_ossos_bias::export()),
        (
            "p9-2017-resonance-hopping",
            p9_2017_resonance_hopping::export(),
        ),
        ("p9-2018-bp519", p9_2018_bp519::export()),
        (
            "p9-2018-chaotic-dynamics",
            p9_2018_chaotic_dynamics::export(),
        ),
        ("p9-2018-kuiper-belt", p9_2018_kuiper_belt::export()),
        ("p9-2018-low-perihelion", p9_2018_low_perihelion::export()),
        ("p9-2018-resonance", p9_2018_resonance::export()),
        (
            "p9-2018-secular-dynamics",
            p9_2018_secular_dynamics::export(),
        ),
        ("p9-2018-wise-search", p9_2018_wise_search::export()),
        ("p9-2019-clustering", p9_2019_clustering::export()),
        (
            "p9-2019-ossos-scattering",
            p9_2019_ossos_scattering::export(),
        ),
        ("p9-2019-review", p9_2019_review::export()),
        ("p9-2019-selfgrav-disk", p9_2019_selfgrav_disk::export()),
        ("p9-2019-tess-yield", p9_2019_tess_yield::export()),
        (
            "p9-2020-clement-precession",
            p9_2020_clement_precession::export(),
        ),
        ("p9-2020-des-isotropy", p9_2020_des_isotropy::export()),
        (
            "p9-2020-resonance-hopping",
            p9_2020_resonance_hopping::export(),
        ),
        (
            "p9-2020-secular-octupole",
            p9_2020_secular_octupole::export(),
        ),
        ("p9-2020-tess-shiftstack", p9_2020_tess_shiftstack::export()),
        ("p9-2021-act-mm", p9_2021_act_mm::export()),
        ("p9-2021-des-catalog", p9_2021_des_catalog::export()),
        (
            "p9-2021-detached-inclinations",
            p9_2021_detached_inclinations::export(),
        ),
        ("p9-2021-napier-critique", p9_2021_napier_critique::export()),
        ("p9-2021-oort-cloud", p9_2021_oort_cloud::export()),
        ("p9-2021-orbit", p9_2021_orbit::export()),
        ("p9-2021-perihelion-gap", p9_2021_perihelion_gap::export()),
        ("p9-2021-stability", p9_2021_stability::export()),
        ("p9-2021-ztf", p9_2021_ztf::export()),
        ("p9-2022-des", p9_2022_des::export()),
        ("p9-2022-iras-candidate", p9_2022_iras_candidate::export()),
        ("p9-2022-uranus-tilt", p9_2022_uranus_tilt::export()),
        ("p9-2023-lsst-strategy", p9_2023_lsst_strategy::export()),
        ("p9-2023-mond-efe", p9_2023_mond_efe::export()),
        (
            "p9-2024-neptune-crossing",
            p9_2024_neptune_crossing::export(),
        ),
        ("p9-2024-oort-selfgrav", p9_2024_oort_selfgrav::export()),
        ("p9-2024-panstarrs", p9_2024_panstarrs::export()),
        (
            "p9-2024-primordial-alignment",
            p9_2024_primordial_alignment::export(),
        ),
        ("p9-2024-siraj-orbit", p9_2024_siraj_orbit::export()),
        (
            "p9-2025-akari-refutation",
            p9_2025_akari_refutation::export(),
        ),
        ("p9-2025-clustering", p9_2025_clustering::export()),
        ("p9-2025-iras-akari", p9_2025_iras_akari::export()),
        ("p9-2025-new-discoveries", p9_2025_new_discoveries::export()),
        ("p9-2025-parallax-search", p9_2025_parallax_search::export()),
        ("p9-2025-perturbation", p9_2025_perturbation::export()),
        ("p9-2025-planet-y", p9_2025_planet_y::export()),
        ("p9-2025-ps1-holman", p9_2025_ps1_holman::export()),
        ("p9-2025-russell-albedo", p9_2025_russell_albedo::export()),
        ("p9-2025-simons-forecast", p9_2025_simons_forecast::export()),
        ("p9-2025-stacking", p9_2025_stacking::export()),
        ("p9-2025-stellar-flybys", p9_2025_stellar_flybys::export()),
        ("p9-2026-alpha-slope", p9_2026_alpha_slope::export()),
        (
            "p9-2026-apsidal-clustering",
            p9_2026_apsidal_clustering::export(),
        ),
        (
            "p9-2026-cluster-inclinations",
            p9_2026_cluster_inclinations::export(),
        ),
        ("p9-2026-flyby-evolution", p9_2026_flyby_evolution::export()),
        (
            "p9-2026-iorio-precession",
            p9_2026_iorio_precession::export(),
        ),
        (
            "p9-2026-stellar-companions",
            p9_2026_stellar_companions::export(),
        ),
    ])
}
