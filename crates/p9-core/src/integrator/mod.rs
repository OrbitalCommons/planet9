pub mod bulirsch_stoer;
pub mod hybrid;
pub mod kepler_step;
pub mod kick;
pub mod whm;

/// Result of a single integration step.
pub struct StepResult {
    /// Number of particles removed this step
    pub removed: usize,
}
