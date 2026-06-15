//! Reproduction of Cowan, Holder & Kaib (2016), "Cosmologists in Search of
//! Planet Nine: the Case for CMB experiments" (arXiv:1602.05963).
//!
//! The paper's headline: a cold (~40 K) Planet Nine's spectral energy
//! distribution peaks in the far-IR/sub-mm, so millimetre CMB experiments and
//! far-IR surveys are well suited to detect it. Relative to representative
//! facility sensitivities the body is far brighter at mm wavelengths than in
//! the optical/near-IR, where it shines only by reflected sunlight.
//!
//! See [`sed`] for the blackbody SED, the Wien peak, the thermal/reflected
//! crossover wavelength, and band-by-band detection margins.

pub mod sed;

pub use sed::P9Sed;
