# Planet Nine — *The Whole Sky* (manim film)

A long-form [manim](https://www.manim.community/) film explaining the physics
behind every paper reproduced in this workspace, assembled from one self-contained
scene module per topic. See issue #185 for the full design.

## Layout

```
animations/
  p9_manim/        shared style + reusable mobjects (theme, layout, orbits, dataio)
  scenes/
    preface/       Act 0 — orbital-mechanics + observational-astronomy primer (P00–P13)
    papers/<crate>/scene.py   one module per reproduced paper
  manifest.yaml    ordered scene list (preface + papers)
  Makefile         render helpers
```

The film is **data-driven**: physics numbers come from the Rust crates (e.g.
`figures/viability.json` via `cargo run -p p9-viability`), loaded by
`p9_manim.dataio`. Python only draws what Rust computed.

## Setup

```bash
pip install --user -r animations/requirements.txt   # manim 0.20.x
# needs system ffmpeg + LaTeX (latex/dvisvgm) + cairo, all standard.
```

## Render

```bash
cd animations
make preface              # render the whole preface at 480p15 (Q=l)
make file  FILE=scenes/preface/preface_a_foundations.py
make scene FILE=scenes/preface/preface_a_foundations.py SCENE=P01Ellipse
make preface Q=h          # 1080p60
```

Rendered video/image artifacts land in `media/` (gitignored).

## Preface curriculum (P00–P13)

Newton's ellipse → Kepler's laws → orbital elements (incl. ϖ) → outer-system
geography → the clustering clue → precession → resonance → how P9 sculpts the
belt → reflected-light & thermal detection → surveys & selection bias →
indirect dynamical fingerprints → the current mass×distance viable region
(bridges to the survey act). Targets a physics-literate audience new to
astro/orbital mechanics.
