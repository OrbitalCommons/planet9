//! Reproduction of Sefilian & Touma (2019)
//! "Shepherding in a Self-gravitating Disk of Trans-Neptunian Objects"
//! (arXiv:1804.06859).
//!
//! Headline claim reproduced here: a massive (~10 Earth-mass), eccentric,
//! apsidally-aligned debris disk in the 100-1000 AU region can secularly
//! SHEPHERD test ETNOs into apsidal confinement (libration of the longitude of
//! perihelion about the disk apse) WITHOUT any planet. The mechanism is
//! Laplace-Lagrange secular forcing by the disk's self-gravity.
//!
//! # What is computed (no hard-coded answers)
//!
//! - The massive disk is modeled as a set of confocal, apsidally-aligned
//!   eccentric rings (`disk`), total mass M_disk.
//! - Orbit-averaged Laplace-Lagrange secular theory (`secular`) gives, for a
//!   test particle at a ~ 250 AU: the FORCED eccentricity e_forced and forced
//!   apse, the disk-induced apsidal precession rate A = dvarpi/dt, and the
//!   secular timescale 2 pi / A. The Laplace coefficients are built in
//!   `laplace` (softened for the disk's finite radial width); the precession
//!   sign/magnitude is cross-checked against `p9_core::analysis::secular`
//!   numerical ring averaging in the tests.
//! - The secular phase portrait (`portrait`) shows a libration island for
//!   M_disk ~ 10 M_earth (Delta varpi confined) that disappears as
//!   M_disk -> 0 (free circulation).
//!
//! # Pinned results (see module tests)
//!
//! - Precession rate scales LINEARLY with M_disk (5 vs 10 M_earth: factor 2).
//! - For M_disk ~ 10 M_earth at a = 250 AU the secular precession period is of
//!   order 1e8-1e9 yr (computed ~1.7e8 yr at the fiducial 0.2 softening; see
//!   the residual note below).
//! - A libration island exists in (e, Delta varpi); whether a given orbit is
//!   PHYSICALLY shepherded (libration completed within the 4.5 Gyr age) grows
//!   with M_disk and vanishes near zero mass, because the precession period
//!   scales as 1/M_disk.
//!
//! # Residuals / honest caveats
//!
//! - **Forced apse sign.** With a uniform-eccentricity, single-signed-density
//!   disk the forced equilibrium comes out APSIDALLY ALIGNED (Delta varpi = 0)
//!   with the disk, not anti-aligned. The paper's anti-aligned confinement
//!   arises from the disk's eccentricity/density gradients (sign of the
//!   cross-term B), which this minimal uniform-e model does not reproduce; the
//!   shepherding-into-a-libration-island mechanism itself is reproduced. The
//!   sign is documented, not asserted as anti-aligned.
//! - **Timescale normalization.** The ~1e8-1e9 yr secular period requires a
//!   finite disk softening (~0.2 of the local a, a dynamically warm disk's
//!   radial/vertical thickness). A razor-thin disk (softening -> 0) precesses
//!   ~10x faster. The softening is a labelled disk-thickness parameter, not a
//!   knob tuned to a single number.

pub mod cross_check;
pub mod disk;
pub mod laplace;
pub mod portrait;
pub mod secular;
