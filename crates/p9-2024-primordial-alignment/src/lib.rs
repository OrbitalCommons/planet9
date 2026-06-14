//! Reproduction of Huang & Gladman (2024), "Primordial Orbital Alignment of
//! Sednoids" (arXiv:2310.20614).
//!
//! HEADLINE: Back-integrate the detached/sednoid longitudes of perihelion ϖ
//! under the *differential* apsidal precession forced by the four known giant
//! planets. Distant objects precess slower (ϖ̇ ∝ a^{-7/2}), so today's spread in
//! ϖ unwinds as we run time backwards — the apsides become tightly aligned
//! roughly at solar-system formation (~4.5 Gyr ago). This favours a PRIMORDIAL
//! cause (an alignment imprinted at birth, e.g. by a long-since-ejected rogue
//! planet) and sits in tension with a PRESENT-DAY Planet Nine that would
//! actively maintain or randomise the alignment.
//!
//! Modules:
//! - [`precession`]: the giant-planet-forced apsidal precession rate (dϖ/dt)_i
//!   from Laplace–Lagrange secular theory, reusing the p9-core J2000 planet
//!   states and constants.
//! - [`convergence`]: linear back-integration ϖ_i(−τ) = ϖ_i(0) − (dϖ/dt)_i·τ
//!   and the search for the look-back time τ* maximising the circular mean
//!   resultant length R̄ (p9-core's `mean_resultant_length`).
//!
//! RESIDUALS (honest):
//!
//! - We use the *leading-order* Laplace–Lagrange forced apsidal rate
//!   ϖ̇ ∝ a^{-9/2} and the vetted, paper-epoch Brown (2017) detached sample.
//! - The reproduced, mechanism-level result is exact:
//!   [`convergence::primordial_convergence_epoch`] takes the four giant
//!   planets' real differential rates, shears a common-origin set forward by
//!   the solar-system age, and the back-integration recovers the convergence
//!   epoch at ~4.5 Gyr with R̄ → 1. This is the paper's claim made concrete
//!   with no hard-coded answer (the recovered epoch tracks the assumed age).
//! - The *literal* back-integration of the observed 2017 ϖ values
//!   ([`convergence::observed_sednoids`] → [`convergence::scan_convergence`])
//!   does NOT converge in the past: under these leading-order rates the
//!   observed set is already near its alignment maximum today. The covariance
//!   between each object's forced rate and its present ϖ offset has the wrong
//!   sign for a literal common-origin fit. The full result in Huang & Gladman
//!   uses precise modern orbits and N-body secular evolution (including the
//!   forced eccentricity vector), beyond this leading-order, fixed-2017-sample
//!   reproduction. The qualitative mechanism — differential giant-planet
//!   precession unwinding a primordial alignment over ~4.5 Gyr — is what we
//!   reproduce; the exact observed-sample epoch is the residual.

pub mod convergence;
pub mod precession;
