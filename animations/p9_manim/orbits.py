"""Reusable orbital-mechanics mobjects: Sun, Kepler ellipses (focus at origin),
apse pointers, points at a given true anomaly, and ETNO swarms.

All sizes are in *scene units* (the caller scales AU -> scene units). Each
ellipse is built centred then shifted so the Sun-focus sits at the origin, and
rotated by ``varpi`` so the perihelion points in the chosen direction.
"""
import numpy as np
from manim import (
    Arrow,
    Dot,
    Ellipse,
    VGroup,
    rotate_vector,
)

from . import theme as T


def sun(radius=0.16):
    """The Sun at the origin."""
    return Dot(point=np.zeros(3), radius=radius, color=T.SUN).set_z_index(5)


def orbit_point(a, e, nu, varpi=0.0):
    """Cartesian scene-position of true anomaly ``nu`` (rad) on orbit (a, e),
    with the Sun at the origin and perihelion rotated by ``varpi``."""
    r = a * (1.0 - e * e) / (1.0 + e * np.cos(nu))
    p = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    return rotate_vector(p, varpi)


def ellipse_orbit(a, e, color=T.BLUE, varpi=0.0, stroke_width=T.ORBIT_STROKE, opacity=1.0):
    """A Kepler ellipse with the Sun-focus at the origin, perihelion along
    direction ``varpi`` (rad from +x)."""
    b = a * np.sqrt(1.0 - e * e)
    c = a * e
    ell = Ellipse(width=2 * a, height=2 * b, color=color, stroke_width=stroke_width)
    ell.set_stroke(opacity=opacity)
    # centre is offset from the focus by c toward aphelion (-x before rotation)
    ell.shift(np.array([-c, 0.0, 0.0]))
    ell.rotate(varpi, about_point=np.zeros(3))
    return ell


def apse_arrow(a, e, varpi=0.0, color=T.ORANGE, scale=1.0):
    """Arrow from the Sun toward perihelion: the longitude-of-perihelion pointer."""
    tip = orbit_point(a, e, 0.0, varpi) * scale
    return Arrow(start=np.zeros(3), end=tip, color=color, buff=0, stroke_width=T.THIN + 1)


def body_dot(a, e, nu, varpi=0.0, color=T.GREEN, radius=0.06):
    return Dot(point=orbit_point(a, e, nu, varpi), radius=radius, color=color)


def etno_swarm(specs, color=T.GREEN, stroke_width=1.6, opacity=0.8):
    """A VGroup of ellipse orbits from a list of (a, e, varpi) tuples."""
    g = VGroup()
    for a, e, varpi in specs:
        g.add(ellipse_orbit(a, e, color=color, varpi=varpi, stroke_width=stroke_width, opacity=opacity))
    return g


def apse_arrows(varpis, length=2.5, color=T.GREEN, stroke_width=3):
    """A VGroup of apsidal-direction arrows from the Sun, one per longitude of
    perihelion in ``varpis`` (radians). Used by the clustering scenes."""
    g = VGroup()
    for v in varpis:
        tip = length * np.array([np.cos(v), np.sin(v), 0.0])
        g.add(Arrow(np.zeros(3), tip, color=color, buff=0, stroke_width=stroke_width))
    return g


def real_etno_swarm(which="kbo", display_a=2.7, color=T.GREEN, stroke_width=1.9, opacity=0.85):
    """A VGroup of the REAL observed objects' orbits (true e and ϖ), plus the
    measured clustering R-bar. ``which`` = 'kbo' or 'etno'. Semi-major axes are
    range-compressed to scene units so every orbit is visible; eccentricity and
    longitude of perihelion are the real measured values. Returns (group, r_bar).
    """
    from . import dataio

    objs, r_bar = dataio.real_etnos(which)
    g = VGroup()
    if not objs:
        return g, r_bar
    amax = max(o["a"] for o in objs)
    for o in objs:
        a_disp = display_a * (o["a"] / amax) ** 0.5
        g.add(ellipse_orbit(a_disp, min(o["e"], 0.86), color=color,
                            varpi=o["varpi_rad"], stroke_width=stroke_width, opacity=opacity))
    return g, r_bar


def mean_resultant_length(varpis):
    """Circular mean resultant length R-bar in [0, 1]; 0 random, 1 aligned."""
    v = np.asarray(varpis)
    return float(np.hypot(np.mean(np.cos(v)), np.mean(np.sin(v))))


def _solve_eccentric_anomaly(e, M):
    """Newton solve of Kepler's equation M = E - e sin E."""
    M = (M + np.pi) % (2 * np.pi) - np.pi
    E = M if e < 0.8 else np.pi
    for _ in range(60):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E -= dE
        if abs(dE) < 1e-12:
            break
    return E


def nu_from_mean_anomaly(e, M):
    """True anomaly (rad) from mean anomaly M (rad) for an elliptic orbit.

    This gives the physically correct, non-uniform angular speed (fast at
    perihelion, slow at aphelion) when M advances linearly in time.
    """
    E = _solve_eccentric_anomaly(e, M)
    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0),
    )
    return nu
