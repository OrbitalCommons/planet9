# Film style guide: one paper, three files

Every reproduced paper appears in the film as three beats:

1. **Claim** (auto-generated) – the date on the publication ribbon, the citation,
   the headline and the claim under test.
2. **Scene** (`scenes/papers/<crate>/scene.py`) – the physics and the paper's
   actual result, drawn from numbers the reproduction crate computed.
3. **Adds** (auto-generated) – what the field knew after the paper that it did
   not know before, shown on the gauge the paper moves, with the workspace's
   reproduced number next to the published one.

Beats 1 and 3 are generated from the ledger entry. A paper therefore owns exactly
three files:

| file | what it holds |
|---|---|
| `animations/ledger/<crate>.yaml` | citation, date, kind, headline, claim, adds, gauge payload, result readout |
| `crates/p9-anim-data/src/papers/<crate_with_underscores>.rs` | `pub fn export() -> Value`: every number the scene and ledger draw |
| `animations/scenes/papers/<crate_with_underscores>/scene.py` | the scene |

`scenes/papers/p9_2021_ztf/` is the reference implementation of all three.

## Ground rules

- **Python draws, Rust computes.** Every plotted series and every quoted number
  comes from the paper's export in `anim.json -> papers -> <crate>`, produced by
  calling the reproduction crate's own functions. Published values appear only
  as labelled comparisons (`published:` in the ledger, or a "paper: X" label).
  No synthetic stand-in data, no `rng.normal` decoration posing as results, and
  no fallback constants for a missing export: if the export is missing the
  scene should fail loudly.
- **Show the result, not a cartoon of it.** Prefer the real sky map, the real
  sample, the computed curve or histogram over a schematic bar or arrow. A
  schematic is fine for explaining a *mechanism*; it must then be followed by
  the computed result.
- **Say what is new.** The scene should make the paper's specific contribution
  visible: what was measured, excluded, predicted or refuted, and by how much.
- **Disagreement is content.** If the reproduction lands away from the published
  number, show both. Do not tune, and do not hide it.
- **Citations come from the ledger**, which was checked against arXiv. Where a
  crate's docs name different authors or a different title, the ledger wins.

## Ledger entry

```yaml
crate: p9-2021-ztf
cite: Brown & Batygin (2021)        # verified against arXiv; do not change
title: ...                          # verified against arXiv; do not change
arxiv: '2110.13117'
date: 2021-10                       # arXiv posting month
kind: search                        # evidence discovery critique dynamics indirect
                                    # search forecast physical alternative review
headline: First big bite out of the search    # <= 45 characters
claim: One sentence, <= 170 characters: what the paper sets out to show.
adds: One or two sentences, <= 260 characters: what the field knew afterwards
  that it did not know before. Concrete and quantitative.
track: excluded                     # the gauge shown in the Adds beat, or none
excluded: {label: ZTF, cumulative: '@cumulative'}
result: {label: predicted orbits ruled out, key: ztf, fmt: pct1, published: 0.564}
```

A value written `'@key'` is read from the paper's Rust export. Use it for every
number the crate can compute.

Gauge payloads (a paper may carry several; `track` picks the one displayed):

| payload | fields | meaning |
|---|---|---|
| `clustering` | `sigma`, `stance` (`pro`/`con`), `label` | one-sided Gaussian significance this paper assigns to the orbital clustering; `con` for analyses that find it consistent with selection bias |
| `orbit` | `mass` (M⊕), `a` (AU), `e`, `i` (deg) | the paper's best-fit perturber |
| `excluded` | `cumulative`, `label` | fraction of the Brown & Batygin (2021) predicted orbits ruled out by direct searches up to and including this paper |
| `sample` | `n`, `what` | number of distant objects the inference rests on after this paper |

Only attach a payload the paper really moves. Most dynamics, indirect, physical
and alternative papers have `track: none`; their Adds beat is the text plus the
result readout.

`result.fmt` is one of `pct0 pct1 pct2 sigma deg au mearth mag` or a Python
format spec such as `'{:.2f} km'`. `tolerance` (default 0.05) sets how close the
reproduction must be to the published value to earn the tick.

## Scene anatomy

```python
class Ztf2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]   # fail loudly if missing
        self.add(paper.scene_header(CRATE))           # headline + citation strip
        ...                                           # 2-4 visual beats
        layout.show_takeaway(self, "One sentence.")   # always last
```

- Length 25–50 s. Two to four beats, each one idea, each with a caption held by
  `timing.hold_to_read`.
- No title card and no hypothesis/conclusion text cards inside the scene: the
  Claim and Adds beats carry those.
- Keep the class name listed in `manifest.yaml`.

### Frame and safe zones

The frame is 14.2 wide by 8.0 high, origin at the centre.

| zone | y range | use |
|---|---|---|
| header | above 3.4 | `paper.scene_header` only |
| stage | −2.9 to 3.1 | all figures, legends, equations |
| caption | below −3.1 | `layout.caption`, `layout.show_takeaway` |

Nothing may overlap anything else. In particular keep legends off the data,
equations off the labels, and axis labels clear of the caption zone (raise the
axes with a negative `shift_down` when a takeaway follows).

### Colour meaning

| colour | meaning |
|---|---|
| `P.BLUE` | Planet Nine itself |
| `P.GREEN` | observed objects and measurements |
| `P.RED` | ruled out, excluded, null |
| `P.ORANGE` | perturbation signatures; the ecliptic |
| `P.PURPLE` | surveys and instruments; the galactic plane |
| `P.TEAL` | forecasts; still viable |
| `P.MUTED` | furniture |

### Helpers

| helper | use |
|---|---|
| `sky.SkyMap` | RA/Dec map, east left: `.reference_curves()`, `.dots()`, `.dec_band()`, `.box()`, `.cells()`, `.polyline()`, `.legend()`, `.p(ra, dec)` |
| `widgets.labeled_axes(..., numbers=True, y_rotate=True)` | axes with tick numbers and labels |
| `widgets.histogram`, `widgets.curve`, `widgets.marker_line` | data on axes |
| `orbits.*` | Sun, Kepler ellipses, apse arrows, real ETNO swarms, precession |
| `layout.label` | every small text (kerning-safe) |
| `layout.caption`, `layout.show_takeaway` | caption zone |
| `layout.explain_equation` | term-by-term equation walk; at most one per scene, placed on the stage with nothing else showing behind it |
| `paper.result_readout` | a boxed number |
| `timing.hold_to_read` | reading-paced holds |

Put any new helper a scene needs in that scene's own file.

## Checking a scene

```bash
cargo run --release -p p9-anim-data            # from the repo root: refresh anim.json
cd animations
python tools/contact_sheet.py scenes/papers/<crate>/scene.py <SceneClass> --frames 12 --out /tmp/<crate>.png
python tools/contact_sheet.py scenes/companions.py Claim_<crate> --frames 4 --out /tmp/<crate>_claim.png
python tools/contact_sheet.py scenes/companions.py Adds_<crate>  --frames 6 --out /tmp/<crate>_adds.png
```

Open each PNG and look at every frame. A scene is done when nothing overlaps or
runs off the frame, every number on screen traces to the export, the text is
readable at 480p, and the three beats tell one story.
