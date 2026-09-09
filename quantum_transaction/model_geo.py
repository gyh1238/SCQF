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
    fig/real_bs.png           the registered base stations over the same
                              ground, read back out of the picture by
                              `bs_sites`

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

  * APs are placed to carry *equal shares* of the ground a UE may stand on,
    not one per lattice cell, and are then snapped onto a base station that
    is really there.  A lattice would spend APs on the hillside, where they
    cover nobody, while overloading the ones that landed on open ground; and
    equalising by centroid alone is not enough, so the share is imposed as a
    capacity constraint (`_balanced_assign`).  Almost every AP lands on a
    registered mast; the few that cannot reach one fall back to allowed
    ground.
  * UEs are drawn uniformly over that same allowed area rather than over the
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

Run `python model_geo.py` for the footprint and how much of it each mask
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
BS_MAP_FILE = "real_bs.png"

# real_bs.png is a screenshot of the registered-base-station map over the same
# ground, from a different provider at a different zoom and crop.  Both are
# north-up, so a scale and an offset per axis relate them:
#
#     campus_px = (sx, sy) * bs_px + (dx, dy)
#
# and the two scales are separate because they came out different: the source
# crops the two axes differently, and one scale left the masts bunched into
# the left half of a map they cover the whole of.
#
# **This fit is by eye and cannot be anything else.**  Four automatic criteria
# were tried and every one of them failed, which is worth recording so that
# nobody spends the afternoon again:
#
#   edge correlation        0.04          noise
#   mutual information      0.016 nats    noise
#   road-mask correlation   flat optimum over sx 0.9-1.15
#   masts on allowed ground 37-44% across every candidate, against 40% for
#                           points thrown at random
#
# The first two fail because campus.png is a 3D render with extruded buildings
# and real_bs.png a flat near-monochrome map: they share almost no pixel
# statistics.  The last is the interesting one -- at no transform do the masts
# prefer ground the AP mask allows, which says the mask is not evidence about
# where a base station goes.  A rooftop in a 3D render is drawn displaced from
# its own footprint, and a good share of real masts stand at the roadside,
# which the mask excludes.
#
# So: change these numbers if the overlay looks wrong to you, and regenerate
# `fig/bs_overlay.png` with `python model_geo.py` to see what you did.  The `y`
# scale and the offsets started from the running track, whose centre reads
# (258, 108) on the source and (470, 210) on the campus, and were then set by
# eye against the overlay -- which is the only instrument this fit has.
BS_TO_CAMPUS = (2.50, 1.30, -335.0, 69.6)   # sx, sy, dx, dy
BS_MERGE_PX = 18.0          # several operators register masts at one site
BS_SNAP_MAX = 0.70          # units: how far an AP may reach for a real mast
# 0.70 is a trade, and both ends of it are measured.  Reaching further puts
# more APs on real masts (68% at 0.45, 77% here, 83% at 1.0) and costs the
# balancing that put them where the demand is (share spread cv 0.23, 0.28,
# 0.36) -- and 0.40 is the unbalanced Lloyd this was built to beat.
BS_CACHE = "bs_sites.npy"   # extraction takes seconds; the result never changes

AP_SNAP_MAX = 0.55          # units: give up on a cell with nothing nearer
AP_MIN_SEP = 0.35           # units: two APs may not share one rooftop
_STRIDE_AP = 3              # candidate pixels are subsampled; APs are few
_STRIDE_UE = 2
_STRIDE_DEMAND = 8          # coarser still: this one only has to find centroids

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


def bs_sites():
    """The registered base-station masts, in campus pixel coordinates.

    The source is a screenshot, so the positions have to be read back off the
    picture.  Every mast is one teardrop pin of a fixed size, which makes this
    template matching rather than blob finding -- pins overlap constantly, and
    a blob of blue is as often four masts as one.  Matching alone still misses
    the ones underneath, so a matched pin is erased from the working image and
    the match is run again, until a pass finds nothing: that lifts the count
    from 78 to 135 on this map.  A pin marks its site with its lower tip, not
    its centre, so the anchor sits `tip` below the match.

    Masts within `BS_MERGE_PX` are then merged.  Operators register per
    carrier and per band, so one rooftop can carry four pins; 135 detections
    are 120 distinct sites.

    Result is cached to `BS_CACHE` -- delete it to re-extract.
    """
    path = os.path.join(_HERE, BS_CACHE)
    if "bs" not in _cache:
        if os.path.exists(path):
            _cache["bs"] = np.load(path)
        else:
            _cache["bs"] = _extract_bs()
            np.save(path, _cache["bs"])
    return _cache["bs"]


def _extract_bs():
    """Read the pins out of `BS_MAP_FILE`; see `bs_sites` for the why."""
    from scipy import ndimage as ndi
    from scipy.signal import fftconvolve
    from scipy.cluster.hierarchy import fcluster, linkage

    img = _load(BS_MAP_FILE)[..., :3].astype(int)
    r, g, b = img[..., 0], img[..., 1], img[..., 2]
    pin = ndi.binary_fill_holes(ndi.binary_closing(
        (b > 120) & (b - r > 50) & (b - g > 25), np.ones((3, 3))))

    lab, _ = ndi.label(pin)                      # one un-overlapped pin, to
    size = np.bincount(lab.ravel())[1:]          # measure the shape from
    alone = int(np.argmin(np.where((size > 900) & (size < 1400),
                                   size, 1 << 30))) + 1
    ys, xs = np.where(lab == alone)
    tpl = pin[max(ys.min() - 2, 0):ys.max() + 3,
              max(xs.min() - 2, 0):xs.max() + 3].astype(float)
    th, tw = tpl.shape
    ker = tpl - tpl.mean()
    norm = float((tpl * ker).sum())
    tip = int(np.where(tpl.any(axis=1))[0].max()) - th // 2

    work = pin.astype(float)
    hits = []
    while True:
        score = fftconvolve(work, ker[::-1, ::-1], mode="same") / norm
        found = 0
        while True:
            i = int(np.argmax(score))
            if score.flat[i] < 0.50:
                break
            y, x = divmod(i, score.shape[1])
            hits.append((x, y + tip))
            found += 1
            score[max(y - 13, 0):y + 14, max(x - 13, 0):x + 14] = -9.0
            y0, y1 = max(y - th // 2, 0), min(y - th // 2 + th, work.shape[0])
            x0, x1 = max(x - tw // 2, 0), min(x - tw // 2 + tw, work.shape[1])
            work[y0:y1, x0:x1] = np.minimum(work[y0:y1, x0:x1],
                                            1 - tpl[:y1 - y0, :x1 - x0])
        if found == 0:
            break

    sx, sy, dx, dy = BS_TO_CAMPUS
    xy = np.array(hits, dtype=float) * np.array([sx, sy]) + np.array([dx, dy])
    group = fcluster(linkage(xy, "single"), t=BS_MERGE_PX, criterion="distance")
    return np.array([xy[group == k].mean(axis=0) for k in np.unique(group)])


def bs_sites_units(win):
    """The masts in model units, and only the ones inside the region."""
    px = bs_sites()
    x, y = win.to_world(px[:, 0], px[:, 1])
    xy = np.stack([x, y], axis=1)
    keep = ((xy[:, 0] >= 0) & (xy[:, 0] <= win.g) &
            (xy[:, 1] >= 0) & (xy[:, 1] <= win.g))
    return xy[keep]


def _balanced_assign(d2, cap):
    """Give every AP the same amount of ground, nearest ground first.

    Plain Lloyd is the obvious way to spread APs over the usable area and it
    does not do what is wanted here: a centroidal tessellation makes each AP
    the centre of its own cell, which on a domain shaped like a campus still
    leaves cells differing sevenfold.  The small cells are the APs that end up
    serving nobody, so the equalisation has to be a constraint, not a hope.

    Each (ground, AP) pair is taken in order of distance and accepted if the
    ground is still free and the AP is not yet full.  The result is an exactly
    balanced assignment that stays close to nearest-AP, which is the same
    greedy any capacity-constrained tessellation starts from.
    """
    n, k = d2.shape
    order = np.argsort(d2, axis=None)
    owner = np.full(n, -1, dtype=int)
    count = np.zeros(k, dtype=int)
    left = n
    for flat in order:
        p, a = flat // k, flat % k
        if owner[p] < 0 and count[a] < cap:
            owner[p] = a
            count[a] += 1
            left -= 1
            if left == 0:
                break
    return owner


def place_aps(rng, g, win, jitter=0.18, iters=8):
    """`g^2` APs spread over the ground that has users, then put on real roofs.

    A lattice laid straight over a campus spends APs on whatever the square
    happens to contain.  On this map that is a wooded hillside and a river of
    road: those APs cover nobody, their zones carry no UE, and the partition
    then divides a region a third of which no user can stand in.  Meanwhile
    the APs that did land on open ground are asked to admit far more UEs than
    `W_a` allows, which is the same imbalance seen from the other end.

    So the APs are placed where the demand is.  Demand is the UE-allowed area
    itself -- UEs are drawn uniformly over it, so equal area is equal expected
    load -- and Lloyd's algorithm from the jittered lattice moves each AP to
    the centroid of the ground closest to it.  What comes out is the flat
    statement of a planner's rule: every AP serves about the same share of the
    campus that anyone actually occupies.  `seed` still varies the deployment,
    through the lattice the iteration starts from.

    Each AP is then snapped to the nearest pixel it may legally occupy, which
    is usually a step of nothing -- most open ground is AP-allowed too -- and
    at most a step onto the roof or lot next door.  `AP_MIN_SEP` keeps two of
    them off the same roof.
    """
    cand = _allowed_xy("ap", win, _STRIDE_AP)
    demand = _allowed_xy("ue", win, _STRIDE_DEMAND)
    if len(cand) == 0:
        raise ValueError(f"no AP-allowed pixel in {win}")
    if len(demand) == 0:
        raise ValueError(f"no UE-allowed pixel in {win}")

    gx, gy = np.meshgrid(np.arange(g), np.arange(g))
    ap = np.stack([gx.ravel(), gy.ravel()], axis=1) + 0.5
    ap = ap + rng.normal(0.0, jitter, ap.shape)

    cap = int(np.ceil(len(demand) / len(ap)))
    for _ in range(iters):
        d2 = ((demand[:, None, :] - ap[None, :, :]) ** 2).sum(axis=2)
        owner = _balanced_assign(d2, cap)
        for k in range(len(ap)):
            mine = owner == k
            if mine.any():
                ap[k] = demand[mine].mean(axis=0)

    # Where the balancing wants an AP, put it on a mast that is really there.
    # The registered sites outnumber the APs about five to one, so asking each
    # AP for the nearest unused one costs the deployment very little and buys
    # it every position from the real map.  Where the screenshot has no mast
    # in reach -- it stops short of the bottom of the region, and it is a
    # screenshot, not a survey -- the allowed-ground snap takes over, and
    # `on_real_mast` says which APs ended up on which.
    sites = bs_sites_units(win)
    free = np.ones(len(sites), bool)

    nearest = np.full(len(cand), np.inf)   # distance to the nearest AP placed
    placed = []                            # so far, carried forward rather
    for p in ap:                           # than recomputed: the candidate
        pick = None                        # list runs to six figures.
        if len(sites):
            d = np.linalg.norm(sites - p, axis=1)
            ok = free & (d <= BS_SNAP_MAX)
            if placed:
                sep = np.linalg.norm(sites[:, None, :] - np.array(placed)[None],
                                     axis=2).min(axis=1)
                ok &= sep >= AP_MIN_SEP
            if ok.any():
                j = int(np.argmin(np.where(ok, d, np.inf)))
                free[j] = False
                pick = sites[j]
        if pick is None:
            d = np.linalg.norm(cand - p, axis=1)
            ok = (d <= AP_SNAP_MAX) & (nearest >= AP_MIN_SEP)
            if not ok.any():
                continue
            pick = cand[np.argmin(np.where(ok, d, np.inf))]
        placed.append(pick)
        nearest = np.minimum(nearest, np.linalg.norm(cand - pick, axis=1))
    return np.array(placed, dtype=float).reshape(-1, 2)


def on_real_mast(ap_xy, win, tol=1e-6):
    """Which of `ap_xy` sit on a registered mast rather than a fallback."""
    sites = bs_sites_units(win)
    if len(sites) == 0:
        return np.zeros(len(ap_xy), bool)
    d = np.linalg.norm(ap_xy[:, None, :] - sites[None], axis=2).min(axis=1)
    return d <= tol


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


def basemap(win, pad=0.0, wash=0.38, desat=0.85):
    """The campus image cropped to `win`, for drawing under a panel.

    Returns `(rgb, extent)` with `extent` in model units and the array in
    `origin="upper"` order, so it goes straight into `imshow`.

    Two knobs, and they do different jobs.  `wash` blends toward white, which
    lowers everything at once -- push it far enough to make a marker stand out
    and the streets go with it.  `desat` blends toward the image's own
    luminance instead, which costs the map nothing structurally: buildings,
    roads and parks keep every edge they had, they simply stop competing for
    hue with the zone colours and the markers drawn over them.  Taking most of
    the reduction out of saturation rather than out of contrast is what lets
    the map stay legible while the model on top of it reads first.
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
    lum = rgb @ np.array([0.2126, 0.7152, 0.0722])   # Rec. 709
    rgb = desat * lum[..., None] + (1.0 - desat) * rgb
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


def write_overlay(path=None, g=5):
    """Draw the extracted masts on the campus, so the fit can be checked.

    The alignment cannot be scored automatically -- see `BS_TO_CAMPUS` -- so
    it has to be checkable by eye, and a claim nobody can check is worse than
    no claim.  This is that check, written next to the rasters it relates.
    """
    from PIL import Image, ImageDraw
    win = window_for(g)
    img = Image.fromarray(_load(BASEMAP_FILE)[..., :3]).convert("RGB")
    d = ImageDraw.Draw(img)
    r0 = win.row1 - win.side_px
    d.rectangle([win.col0, r0, win.col0 + win.side_px, r0 + win.side_px],
                outline=(31, 111, 180), width=5)
    for x, y in bs_sites():
        d.ellipse([x - 7, y - 7, x + 7, y + 7], fill=(214, 39, 40),
                  outline=(255, 255, 255), width=2)
    path = path or os.path.join(FIG_DIR, "bs_overlay.png")
    img.save(path)
    return path


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
    win = window_for(5)
    print(f"  registered masts in the region: {len(bs_sites_units(win))} "
          f"({len(bs_sites_units(win)) / 25:.1f} per AP)")
    print(f"  wrote {write_overlay()}")
