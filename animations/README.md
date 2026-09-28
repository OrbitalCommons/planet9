# Planet Nine — *The Whole Sky* (manim film)

A long-form [manim](https://www.manim.community/) film explaining the physics
behind every paper reproduced in this workspace, assembled from one self-contained
scene module per topic. See issue #185 for the full design.

## Layout

```
animations/
  p9_manim/        shared style + reusable mobjects (theme, layout, orbits, sky,
                   gauges, ledger, contribution, dataio)
  ledger/<crate>.yaml         one dated entry per paper: citation, claim, what it
                              added, and the gauge it moves
  scenes/
    preface/       Act 0 — orbital-mechanics + observational-astronomy primer (P00–P13)
    papers/<crate>/scene.py   one module per reproduced paper
    companions.py  auto-generated Claim_/Adds_/Learn_ beats
  tools/contact_sheet.py      render a scene and tile frames for review
  manifest.yaml    ordered scene list (preface + papers by act and date + finale)
  STYLE.md         how a paper is put on screen
  Makefile         render helpers
```

## What each paper contributes

Every paper plays as three beats: its **claim** (dated on a ribbon of the whole
literature), its **scene**, and what it **added** — shown on one of four gauges
that run through the film (significance of the clustering, best-fit orbit,
fraction of predicted orbits ruled out, size of the sample), with the number
reproduced in this workspace set against the published one. See
[STYLE.md](STYLE.md).

The film is **data-driven**: physics numbers come from the Rust crates (collected
into `animations/data/anim.json` by `cargo run --release -p p9-anim-data`), loaded by
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

## Preface curriculum

Newton's ellipse and Kepler's laws (P00–P02) → orbital elements, incl. ϖ (P03)
→ outer-system geography to scale (P04) → the clustering clue (P05) →
precession (P06) → resonance (P07) → how Planet Nine sculpts the belt: the
orbit-averaged ring, the anti-aligned island, 4 Gyr of herding (P08) →
reflected light and thermal emission (P09–P10) → slow movers and selection bias
(P11) → indirect fingerprints: Cassini ranging and the solar obliquity (P12) →
the predicted orbits carved away survey by survey (P13). Targets a
physics-literate audience new to astro/orbital mechanics. Their numbers come
from `crates/p9-anim-data/src/preface/`, one module per scene file.
