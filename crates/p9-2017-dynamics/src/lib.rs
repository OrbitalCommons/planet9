//! Reproduction of Batygin & Morbidelli (2017)
//! "Dynamical Evolution Induced by Planet Nine"
//!
//! Characterizes dynamical mechanisms by which Planet Nine reproduces observed
//! orbital clustering, facilitates perihelion detachment from Neptune, and excites
//! extreme inclinations in long-period centaurs. Key finding: secular-only models
//! are insufficient; mean-motion resonances play a crucial role in long-term KBO
//! survival. The critical transition from secular to resonance-dominated dynamics
//! occurs at P/P_9 > 0.1-0.15.
//!
//! Fiducial parameters: a_9 = 700 AU, e_9 = 0.6, m_9 = 10 Earth masses.
//! Also tested: a_9 = 600 AU, e_9 = 0.5, and masses 5, 10, 20 Earth masses.

pub mod hamiltonian;
pub mod phase_space;
pub mod plots;
pub mod resonance;
pub mod simulation;
