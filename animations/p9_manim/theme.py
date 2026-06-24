"""Shared visual theme for the Planet Nine film.

Tokyo-Night dark palette so stills drop straight into the repo docs/portal.
Importing this module sets the global manim background colour.
"""
from manim import config, Text

# ---- palette ----
BG = "#1a1b26"      # background
FG = "#c0caf5"      # primary text / axes
MUTED = "#565f89"   # gridlines / secondary
BLUE = "#7aa2f7"    # Planet Nine itself
GREEN = "#9ece6a"   # observed / real data (ETNOs, detections)
RED = "#f7768e"     # exclusions / nulls / ruled-out
ORANGE = "#e0af68"  # perturbation signatures (precession, ranging, obliquity)
PURPLE = "#bb9af7"  # surveys / instruments
TEAL = "#7dcfff"    # forecasts / not-yet-excluded
YELLOW = "#e0af68"
SUN = "#f7c873"     # the Sun

config.background_color = BG

# manim's default Text font resolves (via Pango) to a serif whose glyph advances
# are laid out unevenly -- visible as broken kerning at every size. Noto Sans
# kerns cleanly, so make it the default for every Text in the film.
FONT = "Noto Sans"
Text.set_default(font=FONT)

# Consistent stroke / font sizing.
ORBIT_STROKE = 3.0
THIN = 1.6
TITLE_SIZE = 44
BODY_SIZE = 30
SMALL_SIZE = 24
