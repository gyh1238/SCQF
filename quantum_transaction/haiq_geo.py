"""
Campus geography for the instance generator.
============================================
The synthetic generator drops APs on a jittered grid and UEs uniformly over
the square, which is the right model for the scaling claim -- densities have
to be independent of `g` for `g` to be a pure size knob -- but it says
nothing about whether the protocol survives a real deployment, where the
positions a planner may choose are not a square at all.

This module supplies that geometry from three rasters of one campus, all on
the same 1894x1489 pixel grid:

    fig/ap_allowed_mask.png   where an AP may stand: rooftops and open lots,
                              roads and the excluded hillside removed
    fig/ue_allowed_mask.png   where a UE may stand: open ground only, so a
                              UE never appears inside a building
    fig/campus.png            the basemap the two masks were traced from

Pixels become model units through `PX_PER_UNIT`, fixed so that one unit is
one nominal AP spacing.  An instance of size `g` therefore occupies a `g` x
`g` unit window of the campus, and growing `g` covers more ground at the
same density -- the growth axis is preserved, it is just no longer flat.

What changes in an instance built this way:

  * APs keep their one-per-cell density, but each is snapped to the nearest
    pixel it is allowed to occupy.  A cell whose surroundings are all road
    or hillside yields no AP, so the AP count can fall below `g^2`.
  * UEs are drawn uniformly over the *allowed area* rather than over the
    square, so they cluster along streets and courtyards.

Both effects push coverage overlap around, which is what the partitioner and
the boundary report actually respond to.

Run `python haiq_geo.py` to re-derive `CENTER_PX` -- it scans the campus for
the window carrying both open ground and rooftops at every size in use.
"""

import os
from dataclasses import dataclass

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(_HERE, "fig")

AP_MASK_FILE = "ap_allowed_mask.png"
UE_MASK_FILE = "ue_allowed_mask.png"
BASEMAP_FILE = "campus.png"

PX_PER_UNIT = 150.0         # pixels per nominal AP spacing
CENTER_PX = (875.0, 875.0)  # (col, row) the instance window is centred on

AP_SNAP_MAX = 0.55          # units: give up on a cell with nothing nearer
AP_MIN_SEP = 0.35           # units: two APs may not share one rooftop
_STRIDE_AP = 3              # candidate pixels are subsampled; APs are few
_STRIDE_UE = 2

_cache = {}


def _load(name):
    """Read one raster, once per process."""
    if name not in _cache:
        from PIL import Image
        path = os.path.join(FIG_DIR, name)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} is missing; the campus rasters "
                f"({AP_MASK_FILE}, {UE_MASK_FILE}, {BASEMAP_FILE}) have to sit "
                f"in {FIG_DIR} for geo=True to work")
        _cache[name] = np.asarray(Image.open(path))
    return _cache[name]


def _mask(kind):
    """The allowed-placement mask for `kind`, as a boolean array."""
    return _load(AP_MASK_FILE if kind == "ap" else UE_MASK_FILE) > 0


def shape_px():
    """(height, width) of the campus rasters, in pixels."""
    return _mask("ap").shape


def _to_world(col, row):
    """Pixel (col, row) -> model units, y measured up from the bottom row."""
    h, _ = shape_px()
    return np.asarray(col) / PX_PER_UNIT, (h - 1 - np.asarray(row)) / PX_PER_UNIT


def _to_px(x, y):
    """Model units -> pixel (col, row); the inverse of `_to_world`."""
    h, _ = shape_px()
    return np.asarray(x) * PX_PER_UNIT, (h - 1) - np.asarray(y) * PX_PER_UNIT


@dataclass(frozen=True)
class Window:
    """The square patch of campus one instance is built on, in model units."""
    x0: float
    y0: float
    side: float

    @property
    def x1(self):
        return self.x0 + self.side

    @property
    def y1(self):
        return self.y0 + self.side

    @property
    def lo(self):
        return np.array([self.x0, self.y0])


def window_for(g):
    """The `g` x `g` unit window, centred on `CENTER_PX` and kept in bounds.

    The window grows about a fixed point, so a larger `g` is the same campus
    seen wider rather than a different campus.  A `g` past the raster is
    clamped to it, and the caller is told by the returned side length.
    """
    h, w = shape_px()
    span_x, span_y = (w - 1) / PX_PER_UNIT, (h - 1) / PX_PER_UNIT
    side = float(min(g, span_x, span_y))
    cx, cy = _to_world(*CENTER_PX)
    x0 = float(np.clip(cx - side / 2, 0.0, span_x - side))
    y0 = float(np.clip(cy - side / 2, 0.0, span_y - side))
    return Window(x0, y0, side)


def _allowed_xy(kind, win, stride):
    """World coordinates of every allowed pixel inside `win`, subsampled."""
    key = ("xy", kind, win, stride)
    if key not in _cache:
        c0, r1 = _to_px(win.x0, win.y0)
        c1, r0 = _to_px(win.x1, win.y1)
        h, w = shape_px()
        r0, r1 = max(int(np.floor(r0)), 0), min(int(np.ceil(r1)) + 1, h)
        c0, c1 = max(int(np.floor(c0)), 0), min(int(np.ceil(c1)) + 1, w)
        sub = _mask(kind)[r0:r1:stride, c0:c1:stride]
        rr, cc = np.nonzero(sub)
        x, y = _to_world(c0 + cc * stride, r0 + rr * stride)
        _cache[key] = np.stack([x, y], axis=1)
    return _cache[key]


def place_aps(rng, g, win, jitter=0.18):
    """One AP per grid cell, snapped to the nearest pixel it may occupy.

    The nominal grid and its jitter are the synthetic generator's, so the
    intended density is unchanged; the snap is what the campus imposes.  A
    cell with no allowed pixel within `AP_SNAP_MAX`, or none left once
    `AP_MIN_SEP` is honoured, contributes no AP -- that is a planner finding
    nowhere to mount, and it is allowed to show.
    """
    cell = win.side / g
    cand = _allowed_xy("ap", win, _STRIDE_AP)
    if len(cand) == 0:
        raise ValueError(f"no AP-allowed pixel in {win}")

    gx, gy = np.meshgrid(np.arange(g), np.arange(g))
    nominal = win.lo + (np.stack([gx.ravel(), gy.ravel()], axis=1) + 0.5) * cell
    nominal = nominal + rng.normal(0.0, jitter * cell, nominal.shape)

    placed = []
    for p in nominal:
        d = np.linalg.norm(cand - p, axis=1)
        ok = d <= AP_SNAP_MAX * cell
        if placed:
            sep = np.linalg.norm(cand[:, None, :] - np.array(placed)[None], axis=2)
            ok &= sep.min(axis=1) >= AP_MIN_SEP * cell
        if not ok.any():
            continue
        placed.append(cand[np.argmin(np.where(ok, d, np.inf))])
    return np.array(placed, dtype=float).reshape(-1, 2)


def place_ues(rng, n, win):
    """`n` UEs drawn uniformly over the allowed *area*, not over the square.

    Sampling the allowed pixels and dithering inside the one picked is a
    uniform draw over that area, and it costs one pass instead of the
    unbounded rejection loop a sparse window would need.
    """
    cand = _allowed_xy("ue", win, _STRIDE_UE)
    if len(cand) == 0:
        raise ValueError(f"no UE-allowed pixel in {win}")
    idx = rng.integers(0, len(cand), size=n)
    half = 0.5 * _STRIDE_UE / PX_PER_UNIT
    return cand[idx] + rng.uniform(-half, half, size=(n, 2))


def basemap(win, pad=0.0, wash=0.30):
    """The campus image cropped to `win`, for drawing under a panel.

    Returns `(rgb, extent)` with `extent` in model units and the array in
    `origin="upper"` order, so it goes straight into `imshow`.  It is blended
    `wash` of the way to white -- enough that the zone shading over it stays
    the thing being read, little enough that the buildings and streets the
    placement follows are still identifiable.
    """
    img = _load(BASEMAP_FILE)
    h, w = img.shape[:2]
    c0, r1 = _to_px(win.x0 - pad, win.y0 - pad)
    c1, r0 = _to_px(win.x1 + pad, win.y1 + pad)
    r0, r1 = max(int(np.floor(r0)), 0), min(int(np.ceil(r1)) + 1, h)
    c0, c1 = max(int(np.floor(c0)), 0), min(int(np.ceil(c1)) + 1, w)

    rgb = img[r0:r1, c0:c1, :3].astype(float) / 255.0
    if img.shape[2] == 4:                       # composite onto white first
        a = img[r0:r1, c0:c1, 3:4].astype(float) / 255.0
        rgb = rgb * a + (1.0 - a)
    rgb = wash + (1.0 - wash) * rgb

    xa, ya = _to_world(c0, r1 - 1)
    xb, yb = _to_world(c1 - 1, r0)
    return rgb, (float(xa), float(xb), float(ya), float(yb))


def scan_windows(sizes=(5, 8), step=25):
    """Re-derive `CENTER_PX`: the centre that serves every size in `sizes`.

    A window is scored by the product of its AP-allowed, UE-allowed and
    rooftop-only fractions, taken at its worst size, so the winner is one
    that has open ground to put UEs on *and* buildings to mount APs on at
    both ends of the size range -- not one that is simply all lawn.
    """
    ap, ue = _mask("ap"), _mask("ue")
    bld = ap & ~ue
    h, w = ap.shape

    def integral(m):
        return np.pad(m.astype(np.int64).cumsum(0).cumsum(1), ((1, 0), (1, 0)))

    ints = [integral(m) for m in (ap, ue, bld)]

    def frac(i, r0, c0, s):
        return (i[r0 + s, c0 + s] - i[r0, c0 + s]
                - i[r0 + s, c0] + i[r0, c0]) / (s * s)

    sides = [int(round(g * PX_PER_UNIT)) for g in sizes]
    big = max(sides)
    best = None
    for rc in range(big // 2, h - big // 2, step):
        for cc in range(big // 2, w - big // 2, step):
            worst = [min(frac(i, rc - s // 2, cc - s // 2, s) for s in sides)
                     for i in ints]
            score = float(np.prod(worst))
            if best is None or score > best[0]:
                best = (score, cc, rc, worst)
    score, cc, rc, worst = best
    print(f"CENTER_PX = ({cc:.1f}, {rc:.1f})   score {score:.4f}")
    print(f"  worst-case fractions over g={tuple(sizes)}: "
          f"AP-allowed {worst[0]:.2f}, UE-allowed {worst[1]:.2f}, "
          f"rooftop-only {worst[2]:.2f}")
    return (float(cc), float(rc))


if __name__ == "__main__":
    h, w = shape_px()
    print(f"rasters {w}x{h} px = {(w - 1) / PX_PER_UNIT:.2f} x "
          f"{(h - 1) / PX_PER_UNIT:.2f} units at {PX_PER_UNIT:g} px/unit")
    for g in (3, 5, 8):
        win = window_for(g)
        print(f"  g={g}: window x[{win.x0:.2f},{win.x1:.2f}] "
              f"y[{win.y0:.2f},{win.y1:.2f}] side {win.side:.2f}")
    scan_windows()
