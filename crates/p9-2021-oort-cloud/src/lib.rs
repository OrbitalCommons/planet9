//! Reproduction of Batygin & Brown (2021)
//! "Injection of Inner Oort Cloud Objects Into the Distant Kuiper Belt by Planet Nine"
//!
//! This paper investigates how Planet Nine shepherds inner Oort cloud (IOC) objects
//! into the distant Kuiper belt (a > 250 AU). The key finding is that IOC objects
//! injected by P9 show weaker longitude-of-perihelion confinement (~67%) compared
//! to scattered disk objects (~88%), and preferentially populate the a > 2000 AU
//! region of the distant Kuiper belt.
//!
//! The birth cluster environment (Plummer sphere, M ~ 1200 Msun, r ~ 0.35 pc)
//! generates the IOC population through stellar encounters during the Sun's
//! early cluster membership.

pub mod injection_simulation;
pub mod oort_cloud;
pub mod plots;
pub mod population_comparison;
