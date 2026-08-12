"""
Snapshot figure: one region, partitioned and solved, opened up.
==============================================================
Panels, each carrying one checkable claim rather than an illustration:

 (a) the partition the budget produced -- APs coloured by zone, boundary UEs
     marked, and every zone labelled with (UEs, state qubits, two-qubit gates)
     so the panel itself shows that the partitioner respected the budget;

 (b) per-zone accepted laws -- measured histogram against the exact
     exp(lambda J_z) reference, plotted against utility rank so the panel
     stays readable at any zone size.  Ranks outside F_z carry exactly zero
     mass, which is strict feasibility made visible;

 (c) why the boundary report is a joint list rather than one marginal per
     UE.  The first two panels show what marginals lose: for the most
     dependent pair of boundary UEs in a zone, the product of marginals backs
     a combination the joint rules out entirely.  The third shows what that
     loss costs, by running the whole protocol both ways;

 (d) the order decimation fixed the boundary UEs: one point per UE.

The instance comes from the same generator as the scaling figure.  Its seed
is chosen for legibility -- compact zones, no one-UE zones, boundary links
that can be followed -- and the caption states where its utility ratio falls
among the candidate seeds so that a presentation choice cannot pass for a
quality one.  `preview_seeds.py` renders the candidates side by side.

Usage:  python make_fig_snapshot.py
"""

import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch

from haiq_instance import make_instance, utility_scale
from haiq_partition import partition_aps
from haiq_protocol import run_protocol
from haiq_reference import solve_centralized
from haiq_cost import BUDGET_DEFAULT

PANEL_DIR = "fig/panels"

G = 5
SEED = 7           # chosen for legibility; see pick_seed and preview_seeds.py
BETA = 1.5
K_ACCEPT = 200
SHOT_BUDGET = 10_000
FS_TITLE = 8       # one type scale for every panel, so nothing drifts
FS_LABEL = 7.5
FS_TICK = 6.5
FS_LEGEND = 6.2
COL_RATIO = "#1f6fb4"

N_ZONE_PANELS = 4
N_BND_PANELS = 2
ABLATION_CACHE = "ablation_data.npz"
ABLATION_G = (4, 5, 6, 7)
ABLATION_SEEDS = range(4)


def _tint(color, frac):
    """`color` blended `frac` of the way from white, returned opaque.

    Semi-transparent fills survive PDF and SVG but not every PNG viewer: a
    13%-alpha region is stored as full-saturation RGB with alpha 33, so a
    viewer that flattens or ignores the alpha channel shows the saturated
    colour instead of the tint.  Baking the blend keeps the figure background
    transparent while making the fills render identically everywhere.
    """
    r, g, b = to_rgb(color)
    return (1 - frac + frac * r, 1 - frac + frac * g, 1 - frac + frac * b)


def _count(n):
    """Compact shot count: 7.6e4 reads worse than 76k on a small panel."""
    if n >= 1e6:
        return f"{n/1e6:.1f}M"
    if n >= 1e3:
        return f"{n/1e3:.0f}k"
    return f"{n:.0f}"


def _save(fig, stem, formats=("pdf", "svg", "png"), **kw):
    """Every output is written on a transparent background."""
    for ext in formats:
        fig.savefig(f"{stem}.{ext}", transparent=True, **kw)


def run_seed(seed):
    inst = make_instance(g=G, seed=seed)
    opt, _, ok = solve_centralized(inst)
    if not ok:
        return None
    part = partition_aps(inst)
    res = run_protocol(inst, part, beta=BETA, ubar=utility_scale(inst),
                       k_accept=K_ACCEPT, shot_budget=SHOT_BUDGET,
                       rng=np.random.default_rng(1000 + seed), collect=True)
    return (100.0 * res["utility"] / opt, seed, inst, part, res, opt)


def pick_seed(seed=SEED, seeds=range(12)):
    """Build the chosen instance, and say where its quality sits among its peers.

    The seed is chosen for legibility -- compact zones, no degenerate one-UE
    zones, boundary links that can be followed -- which is a presentation
    choice.  To keep that from becoming a quality selection, the utility ratio
    of every candidate is returned as well, so the caption can state where the
    displayed one falls in that spread.  `preview_seeds.py` renders the
    candidates side by side.
    """
    chosen = run_seed(seed)
    peers = [r[0] for s in seeds if (r := run_seed(s))]
    return chosen, peers


def panel_map(ax, inst, part, active, zcol):
    """(a) the partition, annotated with the resources it produced.

    Zones are drawn as territories -- each point of the plane is shaded by
    the zone owning its nearest AP -- so a zone reads as one contiguous
    region with real borders rather than as a cluster of identical discs.
    What actually couples two zones is a UE whose candidate RBs are owned on
    both sides of a border, so every boundary UE is drawn joined to the APs
    that offer it a candidate; those links are the ones that cross a border.
    """
    # The territory shading has to cover the axes exactly, or the frame shows a
    # white band where the grid stops.  Fix the extent first -- square, so that
    # `aspect("equal")` cannot pad one side -- and build the grid on it.
    lo = np.minimum(inst.ue_xy.min(axis=0), inst.ap_xy.min(axis=0)) - 0.30
    hi = np.maximum(inst.ue_xy.max(axis=0), inst.ap_xy.max(axis=0)) + 0.30
    span = float((hi - lo).max())
    mid = 0.5 * (lo + hi)
    lo, hi = mid - span / 2, mid + span / 2
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], 480),
                         np.linspace(lo[1], hi[1], 480))
    d2 = ((gx[..., None] - inst.ap_xy[:, 0]) ** 2
          + (gy[..., None] - inst.ap_xy[:, 1]) ** 2)
    terr = part.zone_of[np.argmin(d2, axis=2)]           # zone owning each point

    for z in range(len(part.zones)):
        if not part.zones[z]:
            continue
        m = (terr == z).astype(float)
        ax.contourf(gx, gy, m, levels=[0.5, 1.5],
                    colors=[_tint(zcol.get(z, "#cccccc"), 0.30)], zorder=0)
        ax.contour(gx, gy, m, levels=[0.5], colors="white", linewidths=1.0,
                   zorder=1)

    # boundary UEs are joined to the APs that offer them a candidate RB
    bset = set(part.boundary.tolist())
    for i in bset:
        for r in inst.cand[i]:
            a = inst.ap_xy[inst.rb_owner[r]]
            ax.plot([inst.ue_xy[i, 0], a[0]], [inst.ue_xy[i, 1], a[1]],
                    color=_tint("#c1440e", 0.55), lw=0.55, zorder=2)

    for z, aps in enumerate(part.zones):
        if not aps:
            continue
        xy = inst.ap_xy[aps]
        ax.scatter(xy[:, 0], xy[:, 1], s=120, marker="s",
                   color=zcol.get(z, "#dddddd"), edgecolor="black",
                   linewidth=0.7, zorder=4)

    for i in range(inst.n_ue):
        x, y = inst.ue_xy[i]
        if i in bset:
            ax.scatter(x, y, s=28, marker="o", facecolor="white",
                       edgecolor="#c1440e", linewidth=1.3, zorder=5)
        else:
            ax.scatter(x, y, s=11, marker="o",
                       color=_tint(zcol.get(part.ue_zones[i][0], "#999999"), 0.85),
                       zorder=3)

    for r in active:
        zone = r["zone"]
        st = part.stats[zone.idx]
        c = inst.ap_xy[zone.aps].mean(axis=0)
        ax.annotate(f"Z{zone.idx}\n{st['n_ue']}u {st['q_state']}q\n{st['n2q']}",
                    (c[0], c[1] + 0.30), fontsize=5.2, ha="center", va="center",
                    zorder=6, linespacing=0.95,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"))

    ax.set_title(f"(a) partition under a {BUDGET_DEFAULT} two-qubit budget\n"
                 f"{inst.n_ap} APs, {inst.n_ue} UEs "
                 f"$\\rightarrow$ {len(active)} zones, "
                 f"$|\\mathcal{{B}}|$={len(part.boundary)} boundary UEs "
                 f"({100*part.boundary_density:.0f}%)", fontsize=9, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.legend(handles=[
        Patch(facecolor="#bbbbbb", edgecolor="black", label="AP (colour = zone)"),
        plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#c1440e",
                   mew=1.3, ms=6, label="boundary UE"),
        plt.Line2D([], [], marker="o", ls="", color="#999999", ms=4,
                   label="interior UE"),
        plt.Line2D([], [], color=_tint("#c1440e", 0.55), lw=0.9,
                   label="UE-to-candidate-AP link")],
        loc="upper left", bbox_to_anchor=(0.0, -0.015), fontsize=7.5, ncol=2,
        frameon=False, handletextpad=0.4, columnspacing=1.1,
        title="shading: zone territory (nearest AP)",
        title_fontsize=7.5, alignment="left")


def pick_zone_panels(active, n):
    """Zones for the (b) row, chosen to span the range rather than repeat it.

    Taking the n largest feasible sets picks four zones that happen to be the
    same size, which says nothing about how the partition varies.  This walks
    |F_z| downwards, takes each distinct value once, and within a tie prefers
    the AP count that has appeared least -- so the row shows a one-AP zone
    beside a three-AP one, and a wide feasible set beside a narrow one.
    """
    cand = sorted(active, key=lambda r: -len(r["rows"]))
    chosen, seen, ap_seen = [], set(), {}
    for r in cand:
        f = len(r["rows"])
        if f in seen:
            continue
        peers = sorted((q for q in cand if len(q["rows"]) == f),
                       key=lambda q: (ap_seen.get(len(q["zone"].aps), 0),
                                      len(q["zone"].aps)))
        pick = peers[0]
        chosen.append(pick)
        seen.add(f)
        ap_seen[len(pick["zone"].aps)] = ap_seen.get(len(pick["zone"].aps), 0) + 1
        if len(chosen) == n:
            break
    return chosen


def panel_zone_law(ax, r, zcol, first):
    """(b) accepted law of one zone against the exact Gibbs reference."""
    zone = r["zone"]
    rank = np.argsort(-r["uref"])                    # high utility first
    ref = r["pref"][rank]
    key = {tuple(int(x) for x in row): t for t, row in enumerate(r["rows"][rank])}
    meas = np.zeros(len(ref))
    for row, wt in zip(r["codes"], r["weights"]):
        meas[key[tuple(int(x) for x in row)]] += wt
    if meas.sum() > 0:
        meas /= meas.sum()

    xs = np.arange(len(ref))
    ax.bar(xs, meas, width=0.9, color=_tint(zcol.get(zone.idx, "#888"), 0.85),
           label="measured", zorder=2)
    ax.step(np.concatenate([[-0.5], xs + 0.5]), np.concatenate([[ref[0]], ref]),
            where="pre", color="black", lw=1.0, label=r"$\propto e^{\lambda J_z}$",
            zorder=3)
    pad = max(3, len(ref) * 0.09)
    ax.axvspan(len(ref) - 0.5, len(ref) + pad, color=_tint("#c1440e", 0.09), lw=0)
    ax.text(len(ref) + pad * 0.55, ax.get_ylim()[1] * 0.5, "infeasible:\nzero mass",
            fontsize=5.8, color="#c1440e", ha="center", va="center")
    ax.set_xlim(-0.8, len(ref) + pad)
    n_ap = len(zone.aps)
    ax.set_title(f"zone Z{zone.idx}: {n_ap} AP{'' if n_ap == 1 else 's'}, "
                 f"$N_z$={zone.n_ue}, $|\mathcal{{F}}_z|$={len(ref)}",
                 fontsize=FS_TITLE)
    ax.set_xlabel("assignment, ranked by $J_z$", fontsize=FS_LABEL)
    if first:
        ax.set_ylabel("accepted probability", fontsize=FS_LABEL)
    ax.legend(fontsize=FS_LEGEND, framealpha=0.9, loc="upper right")
    ax.tick_params(labelsize=FS_TICK)


def strongest_pair(rep):
    """The most dependent pair of boundary UEs in a zone's report.

    Returns (tv, ka, kb, joint, product) with the empirical joint over that
    pair and the outer product of its own marginals, both taken from the
    retained draws with their reconstruction weights.
    """
    zone, bl = rep["zone"], rep["zone"].boundary_local
    codes, w = rep["codes"], rep["weights"]
    if len(bl) < 2 or len(codes) == 0:
        return None
    best = None
    for a in range(len(bl)):
        for b in range(a + 1, len(bl)):
            ka, kb = bl[a], bl[b]
            j = np.zeros((len(zone.cand[ka]), len(zone.cand[kb])))
            for row, wt in zip(codes, w):
                j[row[ka], row[kb]] += wt
            if j.sum() <= 0:
                continue
            j = j / j.sum()
            prod = np.outer(j.sum(axis=1), j.sum(axis=0))
            tv = 0.5 * float(np.abs(j - prod).sum())
            if best is None or tv > best[0]:
                best = (tv, ka, kb, j, prod)
    return best


def panel_boundary(ax, item, inst, reports, zcol, first):
    """(c) what a joint report carries that a marginal exchange would not.

    The paper's boundary message is a list of retained *joint* draws, and it
    says the product of single-UE marginals is only a one-pass consensus
    estimate: the correlations among boundary UEs live in the joint list, not
    in the product.  The panel puts the two side by side for the most
    dependent pair of boundary UEs in a zone.  The bar to read first is the
    red one -- a combination the joint rules out entirely, which the product
    still backs with a fifth of its belief.
    """
    ri, (tv, ka, kb, joint, prod) = item
    rep = reports[ri]
    zone = rep["zone"]
    ua, ub = zone.ue[ka], zone.ue[kb]
    na, nb = joint.shape

    combos = [(x, y) for x in range(na) for y in range(nb)]
    xs = np.arange(len(combos))
    jv = np.array([joint[x, y] for x, y in combos])
    pv = np.array([prod[x, y] for x, y in combos])

    # a combination absent from F_z is not rare, it is impossible: the two UEs
    # would break an owned RB or AP limit
    rows = rep["rows"]
    bad = np.array([not ((rows[:, ka] == x) & (rows[:, kb] == y)).any()
                    for x, y in combos])

    order = np.argsort(-jv)                    # largest joint bar first
    combos = [combos[i] for i in order]
    jv, pv, bad = jv[order], pv[order], bad[order]

    ax.bar(xs - 0.19, jv, width=0.36, color=_tint(zcol.get(zone.idx, "#888"), 0.95),
           label="joint report")
    ax.bar(xs[~bad] + 0.19, pv[~bad], width=0.36, facecolor="none",
           edgecolor="#333333", hatch="////", lw=0.7, label="product of marginals")
    ax.bar(xs[bad] + 0.19, pv[bad], width=0.36, facecolor="none",
           edgecolor="#c1440e", hatch="////", lw=1.1)

    # Sorting puts the impossible combination last and the key sits in a
    # corner, so both would land on a bar.  Open a clear band above the tallest
    # bar: key on the left of it, label on the right, joined to its bar by a
    # thin leader so the association survives the distance.
    top = max(jv.max(), pv.max())
    ax.set_ylim(0, top * 1.60)
    for i in np.where(bad)[0]:
        ax.annotate("impossible,\nyet backed", (i + 0.19, top * 1.56),
                    ha="center", va="top", fontsize=FS_LEGEND, color="#c1440e",
                    linespacing=1.0)
        ax.plot([i + 0.19, i + 0.19], [pv[i] + top * 0.02, top * 1.30],
                color="#c1440e", lw=0.6, zorder=1)

    ax.set_title(f"zone Z{zone.idx}, boundary UEs {ua} & {ub}", fontsize=FS_TITLE)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"RB{int(zone.cand[ka][x])}\nRB{int(zone.cand[kb][y])}"
                        for x, y in combos], fontsize=FS_TICK - 0.7)
    ax.set_xlabel("joint choice of the two UEs", fontsize=FS_LABEL)
    if first:
        ax.set_ylabel("probability", fontsize=FS_LABEL)
    ax.legend(fontsize=FS_LEGEND, framealpha=0.9, loc="upper left")
    ax.tick_params(labelsize=FS_TICK)


def ablation_data():
    """Utility with the retained joint list, and with marginals only.

    `use_joint=False` treats each zone as though it had reported one marginal
    per boundary UE: nothing is conditioned after a commitment, so the
    correlations never re-enter.  Everything else -- the sampler, the guard,
    the commitment order -- is identical, so the difference is attributable to
    the report format alone.
    """
    if os.path.exists(ABLATION_CACHE):
        return np.load(ABLATION_CACHE)["data"]
    rows = []
    for g in ABLATION_G:
        for sd in ABLATION_SEEDS:
            inst = make_instance(g=g, seed=sd)
            opt, _, ok = solve_centralized(inst)
            if not ok:
                continue
            part = partition_aps(inst)
            ub = utility_scale(inst)
            kw = dict(beta=BETA, ubar=ub, k_accept=K_ACCEPT,
                      shot_budget=SHOT_BUDGET)
            a = run_protocol(inst, part, rng=np.random.default_rng(5),
                             use_joint=True, **kw)
            b = run_protocol(inst, part, rng=np.random.default_rng(5),
                             use_joint=False, **kw)
            rows.append((100 * b["utility"] / opt, 100 * a["utility"] / opt))
    data = np.array(rows)
    np.savez(ABLATION_CACHE, data=data)
    return data


def panel_ablation(ax, data, first=False):
    """(c3) what dropping the re-conditioning costs.

    One point per instance: the utility reached when each zone's marginals
    are taken once and never revisited, against the utility reached when the
    retained joint list is re-conditioned after every commitment.  Points
    above the diagonal are instances where keeping the joint list helped.
    """
    lo_v, hi_v = data[:, 0], data[:, 1]
    lo = min(data.min(), 96.0) - 0.4
    hi = max(data.max(), 100.0) + 0.4

    ax.plot([lo, hi], [lo, hi], "-", color="#999999", lw=1.0, zorder=1)
    ax.annotate("equal", (hi, hi), textcoords="offset points", xytext=(-4, -4),
                ha="right", va="top", fontsize=FS_LEGEND, color="#777777")
    ax.plot(lo_v, hi_v, "o", ms=5.5, color=COL_RATIO, mec="white", mew=0.7,
            zorder=3)

    won = int((hi_v > lo_v).sum())
    gap = float((hi_v - lo_v).mean())
    # every point sits above the line, so the corner below it stays free
    ax.text(0.96, 0.05, f"+{gap:.1f} pp on average\n{won} of {len(data)} above",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FS_LEGEND + 0.6, color=COL_RATIO, linespacing=1.35)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title("cost of dropping the joint list", fontsize=FS_TITLE)
    ax.set_xlabel("one-pass marginals  [%]", fontsize=FS_LABEL)
    ax.set_ylabel("re-conditioned joint list  [%]", fontsize=FS_LABEL)
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=FS_TICK)


def panel_order(ax, res):
    """(d) the order decimation fixed the boundary UEs.

    One point per boundary UE, not per iteration: the x position is where
    that UE fell in the commitment sequence, so the panel says which
    overlaps were resolved early, not how many times anything was repeated.
    """
    conf = [t["conf"] for t in res["trace"]]
    ax.plot(np.arange(1, len(conf) + 1), conf, ".", ms=4, color="#1f6fb4")
    ax.set_ylim(0.4, 1.02)
    ax.set_xlabel("boundary UE, in decimation order", fontsize=FS_LABEL)
    ax.set_ylabel(r"confidence $\max_v b_i(v)$", fontsize=FS_LABEL)
    n_exc = res["exceptions"]
    ax.set_title(f"each boundary UE fixed once, least\n"
                 f"ambiguous first  ·  {n_exc} re-sample"
                 f"{'' if n_exc == 1 else 's'}", fontsize=FS_TITLE)
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=FS_TICK)


def main():
    (ratio, seed, inst, part, res, opt), all_ratios = pick_seed()
    reports = res["reports"]
    active = [r for r in reports if r["zone"].n_ue > 0]
    holders = res["holders"]
    print(f"snapshot: g={G} seed={seed} UEs={inst.n_ue} zones={len(active)} "
          f"ratio={ratio:.2f}% (peers {min(all_ratios):.2f}-{max(all_ratios):.2f}%, "
          f"mean {np.mean(all_ratios):.2f}%)")

    cmap = plt.get_cmap("tab20")
    zcol = {z: cmap(i % 20) for i, z in
            enumerate(sorted(r["zone"].idx for r in active))}

    fig = plt.figure(figsize=(13.6, 6.4))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.95, 1, 1, 1, 1],
                          height_ratios=[1, 1])

    panel_map(fig.add_subplot(gs[:, 0]), inst, part, active, zcol)

    zone_panels = pick_zone_panels(active, N_ZONE_PANELS)
    for j, r in enumerate(zone_panels):
        panel_zone_law(fig.add_subplot(gs[0, 1 + j]), r, zcol, j == 0)

    pairs = [(ri, best) for ri, rep in enumerate(reports)
             if (best := strongest_pair(rep)) is not None]
    pairs.sort(key=lambda it: -it[1][0])
    pairs = pairs[:N_BND_PANELS]
    for j, item in enumerate(pairs):
        panel_boundary(fig.add_subplot(gs[1, 1 + j]), item, inst, reports,
                       zcol, j == 0)
    abl = ablation_data()
    panel_ablation(fig.add_subplot(gs[1, 1 + N_BND_PANELS]), abl)
    panel_order(fig.add_subplot(gs[1, 4]), res)

    # No figure title and no caption block: the caption belongs to the
    # document that places the figure, and the numbers behind this instance
    # are printed to stdout and recorded in fig/README.md.
    fig.subplots_adjust(left=0.028, right=0.988, top=0.945, bottom=0.075,
                        wspace=0.34, hspace=0.42)
    os.makedirs("fig", exist_ok=True)
    _save(fig, "fig/fig_snapshot", formats=("pdf", "png"), dpi=200)
    plt.close(fig)
    print("wrote fig/fig_snapshot.pdf and .png")

    # ---- the same panels again, standalone and separately editable -------
    os.makedirs(PANEL_DIR, exist_ok=True)
    # clear this script's own panels first: names carry zone and UE indices,
    # so a changed instance would otherwise leave stale files behind
    for old in glob.glob(f"{PANEL_DIR}/snapshot_*"):
        os.remove(old)
    specs = [("snapshot_a_partition", (6.4, 6.4),
              lambda ax: panel_map(ax, inst, part, active, zcol))]
    for j, r in enumerate(zone_panels):
        specs.append((f"snapshot_b{j+1}_zone_law_Z{r['zone'].idx}", (4.2, 3.2),
                      lambda ax, r=r: panel_zone_law(ax, r, zcol, True)))
    for j, item in enumerate(pairs):
        zi = reports[item[0]]["zone"]
        ua, ub = zi.ue[item[1][1]], zi.ue[item[1][2]]
        specs.append((f"snapshot_c{j+1}_joint_Z{zi.idx}_UE{ua}_{ub}", (4.0, 3.2),
                      lambda ax, item=item: panel_boundary(
                          ax, item, inst, reports, zcol, True)))
    specs.append((f"snapshot_c{N_BND_PANELS+1}_ablation_joint_vs_marginals",
                  (3.4, 3.2), lambda ax: panel_ablation(ax, abl)))
    specs.append(("snapshot_d_decimation_order", (4.2, 3.2),
                  lambda ax: panel_order(ax, res)))

    for stem, size, draw in specs:
        f, a = plt.subplots(figsize=size)
        draw(a)
        f.tight_layout()
        _save(f, f"{PANEL_DIR}/{stem}", bbox_inches="tight", dpi=200)
        plt.close(f)
    print(f"wrote {len(specs)} standalone panels to {PANEL_DIR}/ (pdf + svg + png)")


if __name__ == "__main__":
    main()
