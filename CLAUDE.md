# Planet 9

## Project Overview

Research and analysis related to the Planet Nine hypothesis.

## Key References

- Michael E. Brown (Caltech): https://mikebrown.caltech.edu/
- ORCID: https://orcid.org/0000-0002-8255-0545

## Code reuse & imports

When shared code lives in `p9-core`, call it directly — do not hide it behind
indirection in the consuming crate:

- **No alias / rename wrapper functions.** Don't wrap a core (or other-crate)
  function in a same-signature local function just to rename it or to "avoid
  naming core" (e.g. `pub fn holman_payne_orbit() -> P9Params { brown_batygin_orbit() }`).
  Call the canonical function directly at every site. A wrapper is only
  acceptable if it genuinely adds behavior (binds crate-specific arguments,
  changes units, etc.) — not if its body is a bare forwarding call.
- **No import-forwarding / re-export routers.** Don't add `pub use other::item;`
  (or a thin module that only re-exports) solely to keep old call paths working
  after code moves into core. Repoint each call site to the item's canonical
  location instead.
- **Import symbols; don't fully-qualify at the call site.** Bring a symbol into
  scope with `use p9_core::module::thing;` and call it unqualified (`thing()`),
  rather than writing `p9_core::module::thing()` inline at the call site.
