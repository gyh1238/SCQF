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

**The region is the whole campus, and `g` sets how densely it is covered.**
There is no scale constant to choose: the model's unit is one nominal AP
spacing, and the spacing is whatever divides the campus into `g` cells, so
`g` APs across always means the largest square the rasters hold.  A `g` of 5
therefore reads as 25 APs over the campus, not as a 5-unit window cut out of
it -- nothing of the map is spent on ground the instance never uses.

The rasters are wider than they are tall, and the region is square by
construction (`g` x `g` cells of equal size), so the square is taken at full
height and centred: the left and right strips, about a fifth of the area,
fall outside every instance.  Reaching them would mean a rectangular region,
which buys 21% more ground at the cost of a panel that no longer fits the
composite's left column.

What changes in an instance built this way:

  * APs keep their one-per-cell density, but each is snapped to the nearest
    pixel it is allowed to occupy.  A cell whose surroundings are all road
    or hillside yields no AP, so the AP count can fall below `g^2`.
  * UEs are drawn uniformly over the *allowed area* rather than over the
    square, so they cluster along streets and courtyards.

Both effects push coverage overlap around, which is what the partitioner and
the boundary report actually respond to.

Everything the model measures -- the coverage radius, the AP jitter, the snap
distances below -- is in AP spacings, so it is unchanged by the scale.  What
the scale does change is how much campus one cell contains: at `g=5` a cell
is 298 px of map, at `g=7` it is 213, so a larger `g` sees the same ground in
finer grain rather than seeing more ground.  That is a deliberate difference
from the scaling figure's growth axis, which holds ground density fixed and
stays on the square for exactly that reason.

Run `python haiq_geo.py` for the footprint and how much of it each mask
allows.
"""

import os
from dataclasses import dataclass

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(_HERE, "fig")

AP_MASK_FILE = "ap_allowed_mask.png"
UE_MASK_FILE = "ue_allowed_mask.png"
BASEMAP_FILE = "campus.png"

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


@dataclass(frozen=True)
class Window:
    """The campus, in model units: `[0, g] x [0, g]` over the raster square.

    `ppu` is derived, not chosen -- it is the pixel width of one AP spacing at
    this `g`, and it is the only place the rasters' pixel grid meets the
    model's units.
    """
    g: int
    ppu: float          # pixels per AP spacing
    col0: float         # raster column of model x = 0
    row1: float         # raster row of model y = 0 (rows count downwards)
    side_px: int

    @property
    def side(self):
        return float(self.g)

    @property
    def x0(self):
        return 0.0

    @property
    def y0(self):
        return 0.0

    @property
    def x1(self):
        return float(self.g)

    @property
    def y1(self):
        return float(self.g)

    @property
    def lo(self):
        return np.zeros(2)

    def to_world(self, col, row):
        """Pixel (col, row) -> model units, y measured upwards."""
        return ((np.asarray(col) - self.col0) / self.ppu,
                (self.row1 - np.asarray(row)) / self.ppu)

    def to_px(self, x, y):
        """Model units -> pixel (col, row); the inverse of `to_world`."""
        return (self.col0 + np.asarray(x) * self.ppu,
                self.row1 - np.asarray(y) * self.ppu)


def window_for(g):
    """The campus as a `g` x `g` region: the largest square the rasters hold.

    The footprint does not depend on `g` -- every instance covers the same
    ground.  What `g` sets is how many AP cells that ground is divided into,
    and therefore how many pixels of map one cell contains.
    """
    h, w = shape_px()
    side_px = min(h, w) - 1
    return Window(g=int(g), ppu=side_px / float(g),
                  col0=(w - 1 - side_px) / 2.0,
                  row1=(h - 1 + side_px) / 2.0,
                  side_px=int(side_px))


def _allowed_xy(kind, win, stride):
    """World coordinates of every allowed pixel inside `win`, subsampled."""
    key = ("xy", kind, win, stride)
    if key not in _cache:
        c0, r1 = win.to_px(win.x0, win.y0)
        c1, r0 = win.to_px(win.x1, win.y1)
        h, w = shape_px()
        r0, r1 = max(int(np.floor(r0)), 0), min(int(np.ceil(r1)) + 1, h)
        c0, c1 = max(int(np.floor(c0)), 0), min(int(np.ceil(c1)) + 1, w)
        sub = _mask(kind)[r0:r1:stride, c0:c1:stride]
        rr, cc = np.nonzero(sub)
        x, y = win.to_world(c0 + cc * stride, r0 + rr * stride)
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
    cand = _allowed_xy("ap", win, _STRIDE_AP)
    if len(cand) == 0:
        raise ValueError(f"no AP-allowed pixel in {win}")

    gx, gy = np.meshgrid(np.arange(g), np.arange(g))
    nominal = np.stack([gx.ravel(), gy.ravel()], axis=1) + 0.5
    nominal = nominal + rng.normal(0.0, jitter, nominal.shape)

    # distance to the nearest AP placed so far, carried forward rather than
    # recomputed: the candidate list runs to six figures at the finer scales,
    # and the all-pairs form allocated a matrix that size once per AP.
    nearest = np.full(len(cand), np.inf)
    placed = []
    for p in nominal:
        d = np.linalg.norm(cand - p, axis=1)
        ok = (d <= AP_SNAP_MAX) & (nearest >= AP_MIN_SEP)
        if not ok.any():
            continue
        pick = cand[np.argmin(np.where(ok, d, np.inf))]
        placed.append(pick)
        nearest = np.minimum(nearest, np.linalg.norm(cand - pick, axis=1))
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
    half = 0.5 * _STRIDE_UE / win.ppu
    return cand[idx] + rng.uniform(-half, half, size=(n, 2))


def basemap(win, pad=0.0, wash=0.22):
    """The campus image cropped to `win`, for drawing under a panel.

    Returns `(rgb, extent)` with `extent` in model units and the array in
    `origin="upper"` order, so it goes straight into `imshow`.  It is blended
    `wash` of the way to white -- enough that the zone outlines over it stay
    the thing being read, little enough that the buildings and streets the
    placement follows are still identifiable.
    """
    img = _load(BASEMAP_FILE)
    h, w = img.shape[:2]
    c0, r1 = win.to_px(win.x0 - pad, win.y0 - pad)
    c1, r0 = win.to_px(win.x1 + pad, win.y1 + pad)
    r0, r1 = max(int(np.floor(r0)), 0), min(int(np.ceil(r1)) + 1, h)
    c0, c1 = max(int(np.floor(c0)), 0), min(int(np.ceil(c1)) + 1, w)

    rgb = img[r0:r1, c0:c1, :3].astype(float) / 255.0
    if img.shape[2] == 4:                       # composite onto white first
        a = img[r0:r1, c0:c1, 3:4].astype(float) / 255.0
        rgb = rgb * a + (1.0 - a)
    rgb = wash + (1.0 - wash) * rgb

    xa, ya = win.to_world(c0, r1 - 1)
    xb, yb = win.to_world(c1 - 1, r0)
    return rgb, (float(xa), float(xb), float(ya), float(yb))


def coverage(g=5):
    """How much of the region each mask allows, and what a cell is worth.

    Reported rather than assumed: the AP count can only fall short of `g^2`
    where the mask leaves a cell nothing to snap to, and the UE draw is
    uniform over the allowed area, so these two fractions are what the
    geometry actually hands the generator.
    """
    win = window_for(g)
    h, w = shape_px()
    r0 = int(round(win.row1)) - win.side_px
    c0 = int(round(win.col0))
    sl = (slice(r0, r0 + win.side_px), slice(c0, c0 + win.side_px))
    ap, ue = _mask("ap")[sl], _mask("ue")[sl]
    return dict(side_px=win.side_px, ppu=win.ppu,
                used_area=win.side_px ** 2 / float(h * w),
                ap_allowed=float(ap.mean()), ue_allowed=float(ue.mean()),
                rooftop_only=float((ap & ~ue).mean()))


if __name__ == "__main__":
    h, w = shape_px()
    print(f"rasters {w} x {h} px")
    for g in (4, 5, 6, 7):
        c = coverage(g)
        print(f"  g={g}: region {c['side_px']} px square "
              f"({c['used_area'] * 100:.0f}% of the rasters), "
              f"one AP spacing = {c['ppu']:.0f} px")
    c = coverage(5)
    print(f"  allowed inside the region: AP {c['ap_allowed']:.2f}, "
          f"UE {c['ue_allowed']:.2f}, rooftop-only {c['rooftop_only']:.2f}")
