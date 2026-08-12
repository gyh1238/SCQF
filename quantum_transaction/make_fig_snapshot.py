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

 (c) boundary coordination -- for the most contested boundary UEs, the
     marginals held by each owning zone, their product b_i, and the value
     decimation committed.  This is what the retained joint list buys and
     what a scalar-preference exchange cannot reproduce;

 (d) the decimation order itself: confidence against commitment step.

The instance is drawn from the same generator as the scaling figure and is
the median-quality run of its size, so the picture is typical, not selected.

Usage:  python make_fig_snapshot.py
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle

from haiq_instance import make_instance, utility_scale
from haiq_partition import partition_aps
from haiq_protocol import run_protocol
from haiq_reference import solve_centralized
from haiq_cost import BUDGET_DEFAULT

G = 5
BETA = 1.5
K_ACCEPT = 2000
N_ZONE_PANELS = 4
N_BND_PANELS = 3


def pick_median_seed(seeds=range(8)):
    """Choose the instance whose utility ratio is the median of its size."""
    runs = []
    for s in seeds:
        inst = make_instance(g=G, seed=s)
        opt, _, ok = solve_centralized(inst)
        if not ok:
            continue
        part = partition_aps(inst)
        res = run_protocol(inst, part, beta=BETA, ubar=utility_scale(inst),
                           k_accept=K_ACCEPT, rng=np.random.default_rng(1000 + s),
                           collect=True)
        runs.append((100.0 * res["utility"] / opt, s, inst, part, res, opt))
    runs.sort(key=lambda t: t[0])
    return runs[len(runs) // 2], [r[0] for r in runs]


def panel_map(ax, inst, part, active, zcol):
    """(a) the partition, annotated with the resources it produced."""
    for z, aps in enumerate(part.zones):
        if not aps:
            continue
        col = zcol.get(z, "#dddddd")
        for a in aps:                       # soft blobs merge adjacent APs
            ax.add_patch(Circle(inst.ap_xy[a], 0.46, facecolor=col, alpha=0.22,
                                lw=0, zorder=1))
        xy = inst.ap_xy[aps]
        ax.scatter(xy[:, 0], xy[:, 1], s=120, marker="s", color=col,
                   edgecolor="black", linewidth=0.6, zorder=3)

    bset = set(part.boundary.tolist())
    for i in range(inst.n_ue):
        x, y = inst.ue_xy[i]
        if i in bset:
            ax.scatter(x, y, s=26, marker="o", facecolor="white",
                       edgecolor="#c1440e", linewidth=1.3, zorder=4)
        else:
            ax.scatter(x, y, s=11, marker="o",
                       color=zcol.get(part.ue_zones[i][0], "#999999"),
                       alpha=0.8, zorder=2)

    for r in active:
        zone = r["zone"]
        st = part.stats[zone.idx]
        c = inst.ap_xy[zone.aps].mean(axis=0)
        ax.annotate(f"Z{zone.idx}\n{st['n_ue']}u {st['q_state']}q\n{st['n2q']}",
                    (c[0], c[1] + 0.30), fontsize=5.2, ha="center", va="center",
                    zorder=6, linespacing=0.95,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                              alpha=0.80))

    ax.set_title(f"(a) partition under a {BUDGET_DEFAULT} two-qubit budget\n"
                 f"{inst.n_ap} APs, {inst.n_ue} UEs "
                 f"$\\rightarrow$ {len(active)} zones, "
                 f"$|\\mathcal{{B}}|$={len(part.boundary)} boundary UEs "
                 f"({100*part.boundary_density:.0f}%)", fontsize=9, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    ax.margins(0.03)
    ax.legend(handles=[
        Patch(facecolor="#bbbbbb", edgecolor="black", label="AP (colour = zone)"),
        plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#c1440e",
                   mew=1.3, ms=6, label="boundary UE"),
        plt.Line2D([], [], marker="o", ls="", color="#999999", ms=4,
                   label="interior UE")],
        loc="upper left", bbox_to_anchor=(0.0, -0.015), fontsize=7.5, ncol=3,
        frameon=False, handletextpad=0.4, columnspacing=1.1,
        title="zone labels: UEs, state qubits, two-qubit gates",
        title_fontsize=7.5, alignment="left")


def panel_zone_law(ax, r, zcol, first):
    """(b) accepted law of one zone against the exact Gibbs reference."""
    zone = r["zone"]
    rank = np.argsort(-r["uref"])                    # high utility first
    ref = r["pref"][rank]
    key = {tuple(int(x) for x in row): t for t, row in enumerate(r["rows"][rank])}
    meas = np.zeros(len(ref))
    for row in r["codes"]:
        meas[key[tuple(int(x) for x in row)]] += 1
    meas /= max(len(r["codes"]), 1)

    xs = np.arange(len(ref))
    ax.bar(xs, meas, width=0.9, color=zcol.get(zone.idx, "#888"), alpha=0.85,
           label="measured", zorder=2)
    ax.step(np.concatenate([[-0.5], xs + 0.5]), np.concatenate([[ref[0]], ref]),
            where="pre", color="black", lw=1.0, label=r"$\propto e^{\lambda J_z}$",
            zorder=3)
    pad = max(3, len(ref) * 0.09)
    ax.axvspan(len(ref) - 0.5, len(ref) + pad, color="#c1440e", alpha=0.09, lw=0)
    ax.text(len(ref) + pad * 0.55, ax.get_ylim()[1] * 0.5, "infeasible:\nzero mass",
            fontsize=5.8, color="#c1440e", ha="center", va="center")
    ax.set_xlim(-0.8, len(ref) + pad)
    ax.set_title(f"zone Z{zone.idx}: $N_z$={zone.n_ue}, "
                 f"$|\\mathcal{{F}}_z|$={len(ref)}\n"
                 f"$\\mu_z$={r['mu']:.3f}, $k$={r['rounds']} rounds", fontsize=7.5)
    ax.set_xlabel("assignment, ranked by $J_z$", fontsize=7)
    if first:
        ax.set_ylabel("accepted probability", fontsize=8)
        ax.legend(fontsize=6.2, framealpha=0.9, loc="upper right")
    ax.tick_params(labelsize=6.5)


def panel_boundary(ax, t, inst, reports, holders, zcol, first):
    """(c) competing zone marginals, their product, and the committed value."""
    i = t["ue"]
    n_val = len(inst.cand[i])
    hs = holders[i]
    width = 0.8 / (len(hs) + 1)
    xs = np.arange(n_val)
    for hj, (ri, k) in enumerate(hs):
        rep = reports[ri]
        c = np.bincount(rep["codes"][:, k], minlength=n_val).astype(float)
        c = (c + 0.5) / (c.sum() + 0.5 * n_val)
        ax.bar(xs + hj * width - 0.4 + width / 2, c, width=width * 0.9,
               color=zcol.get(rep["zone"].idx, "#888"), alpha=0.9,
               label=f"$\\pi_{{Z{rep['zone'].idx}}}$")
    ax.bar(xs + len(hs) * width - 0.4 + width / 2, t["belief"],
           width=width * 0.9, color="#222222", alpha=0.9, label="$b_i$")
    ax.axvline(t["value"], color="#c1440e", ls="--", lw=1.3, zorder=0)
    ax.annotate("committed", (t["value"], 0.985), xycoords=("data", "axes fraction"),
                fontsize=6.2, color="#c1440e", ha="center", va="top",
                bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none",
                          alpha=0.85))
    owners = " + ".join(f"Z{reports[ri]['zone'].idx}" for ri, _ in hs)
    ax.set_title(f"boundary UE {i}:  {owners}\nconfidence {t['conf']:.2f}",
                 fontsize=7.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"RB{int(v)}" for v in inst.cand[i]], fontsize=6.5)
    ax.set_xlabel("candidate", fontsize=7)
    if first:
        ax.set_ylabel("marginal / belief", fontsize=8)
    ax.legend(fontsize=5.8, framealpha=0.9, ncol=1, loc="upper left")
    ax.tick_params(labelsize=6.5)


def panel_order(ax, res):
    """(d) decimation commits the least ambiguous boundary UE first."""
    conf = [t["conf"] for t in res["trace"]]
    ax.plot(np.arange(1, len(conf) + 1), conf, ".", ms=4, color="#1f6fb4")
    ax.set_ylim(0.4, 1.02)
    ax.set_xlabel("commitment step", fontsize=7)
    ax.set_ylabel(r"confidence $\max_v b_i(v)$", fontsize=7.5)
    ax.set_title(f"decimation order\n{res['exceptions']} exception re-samples",
                 fontsize=7.5)
    ax.grid(alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=6.5)


def main():
    (ratio, seed, inst, part, res, opt), all_ratios = pick_median_seed()
    reports = res["reports"]
    active = [r for r in reports if r["zone"].n_ue > 0]
    holders = res["holders"]
    print(f"snapshot: g={G} seed={seed} UEs={inst.n_ue} zones={len(active)} "
          f"ratio={ratio:.2f}% (median of {len(all_ratios)} seeds)")

    cmap = plt.get_cmap("tab20")
    zcol = {z: cmap(i % 20) for i, z in
            enumerate(sorted(r["zone"].idx for r in active))}

    fig = plt.figure(figsize=(13.6, 6.9))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.95, 1, 1, 1, 1],
                          height_ratios=[1, 1])

    panel_map(fig.add_subplot(gs[:, 0]), inst, part, active, zcol)

    order = np.argsort([-len(r["rows"]) for r in active])
    for j, r in enumerate([active[i] for i in order[:N_ZONE_PANELS]]):
        panel_zone_law(fig.add_subplot(gs[0, 1 + j]), r, zcol, j == 0)

    contested = sorted([t for t in res["trace"] if len(holders[t["ue"]]) > 1],
                       key=lambda t: t["conf"])[:N_BND_PANELS]
    for j, t in enumerate(contested):
        panel_boundary(fig.add_subplot(gs[1, 1 + j]), t, inst, reports, holders,
                       zcol, j == 0)
    panel_order(fig.add_subplot(gs[1, 4]), res)

    fig.suptitle("A region partitioned to the hardware budget, solved by "
                 "zone-local sampling and classical boundary coordination",
                 fontsize=11.5, y=0.982)
    fig.text(0.008, 0.028,
             f"g={G}, seed={seed} (median of {len(all_ratios)} seeds by utility "
             f"ratio); beta={BETA}, K_z={K_ACCEPT}. Utility {res['utility']:.1f} "
             f"of the centralized strict optimum {opt:.1f} = {ratio:.1f}%; the "
             f"assignment is strictly feasible ({res['feasible']}).",
             fontsize=7, color="#555555")
    fig.text(0.008, 0.008,
             "Zone laws are exact: the accepted branch of the circuit matches "
             "exp(lambda J_z) on F_z to statevector precision (max TVD 2.7e-15 "
             "over the zones checked in haiq_certify.py), so the sampler used "
             "here is the same law the circuit produces.",
             fontsize=7, color="#555555")
    fig.subplots_adjust(left=0.028, right=0.988, top=0.90, bottom=0.115,
                        wspace=0.34, hspace=0.52)
    os.makedirs("fig", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"fig/fig_snapshot.{ext}", dpi=200)
    print("wrote fig/fig_snapshot.pdf and .png")


if __name__ == "__main__":
    main()
