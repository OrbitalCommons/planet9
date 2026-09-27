"""An equatorial sky map (plate carree) with the astronomer's orientation:
right ascension increases to the LEFT, declination up.

    sky = SkyMap()                       # frame, grid, RA/Dec labels
    sky.add_reference_curves()           # ecliptic + galactic plane from Rust
    dots = sky.dots(samples, color=...)  # (ra, dec) points -> Dots
    band = sky.dec_band(-31, 90, ...)    # a declination-limited footprint
    box = sky.box(ra0, ra1, dec0, dec1)  # an RA/Dec rectangle (may wrap RA 0)

Everything is positioned through :meth:`SkyMap.p`, so a map can be scaled or
moved as a group and its helpers stay registered to it.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    DashedVMobject,
    Dot,
    Line,
    Polygon,
    Rectangle,
    VGroup,
    VMobject,
)

from . import dataio, layout
from . import theme as T


class SkyMap(VGroup):
    def __init__(self, width=11.0, dec_range=(-90.0, 90.0), centre=(0.0, -0.35, 0.0),
                 ra_centre=180.0, grid=True, **kwargs):
        super().__init__(**kwargs)
        self.map_w = width
        self.dec_lo, self.dec_hi = dec_range
        self.map_h = width * (self.dec_hi - self.dec_lo) / 360.0
        self.ra_centre = ra_centre
        self.origin = np.array(centre, dtype=float)

        self.frame = Rectangle(width=self.map_w, height=self.map_h, color=T.MUTED,
                               stroke_width=1.4)
        self.frame.set_fill("#16171f", opacity=1.0).move_to(self.origin)
        self.add(self.frame)

        self.graticule = VGroup()
        self.labels = VGroup()
        for ra_h in range(0, 25, 4):
            ra = 15.0 * ra_h
            x = self._x(ra if ra_h < 24 else 359.999)
            if grid and 0 < ra_h < 24:
                self.graticule.add(Line([x, self._y(self.dec_lo), 0], [x, self._y(self.dec_hi), 0],
                                        color=T.MUTED, stroke_width=0.6).set_stroke(opacity=0.35))
            lab = layout.label(f"{ra_h % 24}h", font_size=12, color=T.MUTED)
            lab.move_to([x, self._y(self.dec_lo) - 0.2, 0])
            self.labels.add(lab)
        step = 30
        first = int(np.ceil(self.dec_lo / step) * step)
        for dec in range(first, int(self.dec_hi) + 1, step):
            y = self._y(dec)
            if grid and self.dec_lo < dec < self.dec_hi:
                self.graticule.add(Line([self._x(359.999), y, 0], [self._x(0.0), y, 0],
                                        color=T.MUTED, stroke_width=0.6).set_stroke(opacity=0.35))
            lab = layout.label(f"{dec:+d}°", font_size=12, color=T.MUTED)
            lab.move_to([self.origin[0] - self.map_w / 2 - 0.35, y, 0])
            self.labels.add(lab)
        axis = layout.label("right ascension  (east ←)", font_size=12, color=T.MUTED)
        axis.move_to([self.origin[0], self._y(self.dec_lo) - 0.48, 0])
        self.labels.add(axis)
        self.add(self.graticule, self.labels)
        self._anchor = Dot(self.origin, radius=1e-3).set_opacity(0)
        self._corner = Dot(self.origin + np.array([self.map_w / 2, self.map_h / 2, 0]),
                           radius=1e-3).set_opacity(0)
        self.add(self._anchor, self._corner)

    # ---- coordinate mapping (unit square, then placed by the anchors) -------
    def _u(self, ra):
        """RA (deg) -> horizontal fraction in [-0.5, 0.5], east to the left."""
        d = (float(ra) - self.ra_centre + 180.0) % 360.0 - 180.0
        return -d / 360.0

    def _v(self, dec):
        mid = 0.5 * (self.dec_lo + self.dec_hi)
        return (float(dec) - mid) / (self.dec_hi - self.dec_lo)

    def _x(self, ra):
        return self.origin[0] + self._u(ra) * self.map_w

    def _y(self, dec):
        return self.origin[1] + self._v(dec) * self.map_h

    def p(self, ra, dec):
        """Scene point of (RA, Dec) in degrees, tracking any move/scale of the map."""
        c = self._anchor.get_center()
        half = self._corner.get_center() - c
        return c + np.array([2 * self._u(ra) * half[0], 2 * self._v(dec) * half[1], 0.0])

    def inside(self, dec):
        return self.dec_lo <= dec <= self.dec_hi

    # ---- drawables ----------------------------------------------------------
    def polyline(self, radec, color=T.FG, stroke_width=1.6, opacity=1.0, dashed=False):
        """A curve through (RA, Dec) pairs, broken where it wraps in RA or leaves
        the declination range."""
        runs, cur, last_u = [], [], None
        for ra, dec in radec:
            u = self._u(ra)
            if not self.inside(dec) or (last_u is not None and abs(u - last_u) > 0.5):
                if len(cur) > 1:
                    runs.append(cur)
                cur = []
            if self.inside(dec):
                cur.append(self.p(ra, dec))
            last_u = u
        if len(cur) > 1:
            runs.append(cur)
        g = VGroup()
        for pts in runs:
            m = VMobject(color=color, stroke_width=stroke_width)
            m.set_points_as_corners(pts)
            m.set_stroke(opacity=opacity)
            g.add(DashedVMobject(m, num_dashes=max(8, len(pts) // 3)) if dashed else m)
        return g

    def reference_curves(self):
        """(ecliptic, galactic plane, labels) from the Rust export."""
        s = dataio.section("sky") or {}
        ecl = self.polyline(s.get("ecliptic", []), color=T.ORANGE, stroke_width=1.4,
                            opacity=0.8, dashed=True)
        gal = VGroup(
            self.polyline(s.get("galactic_plane", []), color=T.PURPLE, stroke_width=1.6,
                          opacity=0.7),
            self.polyline(s.get("galactic_b_plus10", []), color=T.PURPLE, stroke_width=0.8,
                          opacity=0.35),
            self.polyline(s.get("galactic_b_minus10", []), color=T.PURPLE, stroke_width=0.8,
                          opacity=0.35),
        )
        return ecl, gal

    def legend(self, items, font_size=13):
        """A row of (label, colour) keys under the map."""
        row = VGroup()
        for text, col in items:
            swatch = Line(LEFT * 0.18, LEFT * -0.18, color=col, stroke_width=3)
            lab = layout.label(text, font_size=font_size, color=T.FG)
            row.add(VGroup(swatch, lab).arrange(buff=0.1))
        row.arrange(buff=0.45)
        row.next_to(self.frame, DOWN, buff=0.62)
        return row

    def dots(self, samples, color=T.TEAL, radius=0.025, opacity=0.8):
        """Dots for samples given as dicts with ra_deg/dec_deg or (ra, dec) pairs."""
        g = VGroup()
        for s in samples:
            ra, dec = (s["ra_deg"], s["dec_deg"]) if isinstance(s, dict) else s
            if self.inside(dec):
                g.add(Dot(self.p(ra, dec), radius=radius, color=color).set_opacity(opacity))
        return g

    def dec_band(self, dec_lo, dec_hi, color=T.PURPLE, opacity=0.13, stroke_width=1.5):
        """The sky between two declinations (a ground survey's reach)."""
        lo, hi = max(dec_lo, self.dec_lo), min(dec_hi, self.dec_hi)
        a, b = self.p(359.999 + self.ra_centre - 180.0, lo), self.p(self.ra_centre - 179.999, hi)
        r = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=color,
                    stroke_width=stroke_width)
        return r.set_fill(color, opacity=opacity)

    def box(self, ra_lo, ra_hi, dec_lo, dec_hi, color=T.PURPLE, opacity=0.13, stroke_width=1.5):
        """An RA/Dec rectangle from ``ra_lo`` eastward to ``ra_hi`` (wrapping
        through RA 0 when ra_hi < ra_lo); returns a VGroup of 1-2 polygons."""
        spans = [(ra_lo, ra_hi)] if ra_hi >= ra_lo else [(ra_lo, 359.999), (0.0, ra_hi)]
        edge = (self.ra_centre + 180.0) % 360.0
        pieces = []
        for lo, hi in spans:
            if lo < edge < hi:
                pieces += [(lo, edge - 1e-3), (edge + 1e-3, hi)]
            else:
                pieces.append((lo, hi))
        g = VGroup()
        for lo, hi in pieces:
            a, b = self.p(lo, max(dec_lo, self.dec_lo)), self.p(hi, min(dec_hi, self.dec_hi))
            poly = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=color,
                           stroke_width=stroke_width)
            g.add(poly.set_fill(color, opacity=opacity))
        return g

    def cells(self, ra_centres, dec_centres, values, color=T.ORANGE, vmax=None, max_opacity=0.9,
              floor=0.02):
        """A gridded map (row-major, dec rows x ra columns) as translucent cells
        whose opacity scales with value / vmax."""
        vals = np.asarray(values, dtype=float).reshape(len(dec_centres), len(ra_centres))
        vmax = float(vmax if vmax is not None else np.nanmax(vals))
        dra = abs(ra_centres[1] - ra_centres[0])
        ddec = abs(dec_centres[1] - dec_centres[0])
        g = VGroup()
        for iy, dec in enumerate(dec_centres):
            if not self.inside(dec):
                continue
            for ix, ra in enumerate(ra_centres):
                f = vals[iy, ix] / vmax if vmax > 0 else 0.0
                if not np.isfinite(f) or f < floor:
                    continue
                a = self.p(ra - dra / 2, dec - ddec / 2)
                b = self.p(ra + dra / 2, dec + ddec / 2)
                cell = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0)
                g.add(cell.set_fill(color, opacity=max_opacity * min(f, 1.0)))
        return g
