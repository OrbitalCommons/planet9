//! Where should a small space telescope image to find Planet Nine, or the
//! distant objects that would stand in for it?
//!
//! The answer is built from four ingredients the workspace already has and
//! one it did not:
//!
//! 1. **Where it could be** ([`prior`]): draws of the predicted orbit, placed
//!    on tonight's sky with their brightness.
//! 2. **What has been done** ([`prior`]): each draw scored by the reproduced
//!    ZTF, DES and Pan-STARRS1 searches.
//! 3. **What will be done anyway** ([`prior`], [`field`]): Rubin's reach,
//!    conceded rather than duplicated.
//! 4. **What this telescope can do** ([`instrument`]): depth against
//!    integration for a read-noise-limited small aperture, linking of a slow
//!    mover, and the cost of a tile.
//! 5. **The Galaxy** ([`crowding`]): every telescope's completeness falls
//!    toward the plane at a rate set by its image quality. That single
//!    consideration decides the answer: the ground surveys are blind where
//!    the predicted orbit crosses the Milky Way, and a sharp space telescope
//!    is not.
//!
//! [`optimize`] then spends a time budget tile by tile and tier by tier in
//! order of probability captured per hour, [`season`] says when each field is
//! at opposition, and [`tno`] counts the distant objects that come for free.

pub mod crowding;
pub mod field;
pub mod instrument;
pub mod optimize;
pub mod prior;
pub mod report;
pub mod season;
pub mod svg;
pub mod tiles;
pub mod tno;
pub mod zones;
