"""Shared visual theme for the Planet Nine film.

Tokyo-Night dark palette so stills drop straight into the repo docs/portal.
Importing this module sets the global manim background colour.
"""
from manim import config

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

# Consistent stroke / font sizing.
ORBIT_STROKE = 3.0
THIN = 1.6
TITLE_SIZE = 44
BODY_SIZE = 30
SMALL_SIZE = 24
