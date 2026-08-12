"""
Scaling figure: what stays bounded as the region grows.
=======================================================
Growth axis: the region is enlarged at fixed AP spacing, UE density,
coverage radius and candidate degree, so the *global* problem grows while
local density does not.  The partitioner then emits more zones of roughly
constant size rather than larger zones.

(a) two-qubit gates of one oracle pass -- the single centralized circuit
    against the largest zone circuit, with the budget the partitioner
    respects;
(b) classical coordination traffic -- the whole region against one zone;
(c) utility of the distributed result as a fraction of the centralized
    strict optimum (MILP).

Read together: the work one processor has to do, and the traffic one zone
has to send, are both unchanged as the region grows, while the centralized
circuit and the total traffic grow with it -- at no cost in solution
quality.  No competing protocol is needed for this; the centralized optimum
enters only as the denominator of (c).

Usage:  python make_fig_scaling.py [--recollect]
"""

import argparse
import glob
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

from haiq_instance import make_instance, utility_scale
from haiq_partition import partition_aps
from haiq_cost import centralized_2q_cost, BUDGET_DEFAULT
from haiq_protocol import run_protocol
from haiq_reference import solve_centralized

CACHE = "scaling_data.npz"
G_VALUES = (3, 4, 5, 6, 7, 8, 9)
SEEDS = tuple(range(6))
BETA = 1.5
K_ACCEPT = 200
SHOT_BUDGET = 10_000


def collect(verbose=True):
    rec = []
    skipped = 0
    for g in G_VALUES:
        for s in SEEDS:
            t0 = time.time()
            inst = make_instance(g=g, seed=s)
            opt, _, ok = solve_centralized(inst)
            if not ok:
                skipped += 1
                if verbose:
                    print(f"  g={g} seed={s}: infeasible instance, skipped")
                continue
            part = partition_aps(inst)
            res = run_protocol(inst, part, beta=BETA, ubar=utility_scale(inst),
                               k_accept=K_ACCEPT, shot_budget=SHOT_BUDGET,
                               rng=np.random.default_rng(1000 + s))
            cen = centralized_2q_cost(inst)
            rec.append((g, s, inst.n_ue, res["n_active_zones"],
                        100.0 * res["utility"] / opt,
                        max(st["n2q"] for st in part.stats),
                        cen["n2q"],
                        max(st["n_qubits"] for st in part.stats),
                        cen["n_qubits"],
                        100.0 * part.boundary_density,
                        res["comm_bits"], res["exceptions"],
                        float(res["feasible"]), res["shots_max"],
                        res["attempts"], res["n_backed_off"], res["ess_min"],
                        res["beta_min"]))
            if verbose:
                print(f"  g={g} seed={s}: UEs={inst.n_ue:3d} zones={res['n_active_zones']:2d} "
                      f"ratio={rec[-1][4]:.2f}% feas={res['feasible']} "
                      f"({time.time()-t0:.1f}s)")
    arr = np.array(rec, dtype=float)
    np.savez(CACHE, data=arr, skipped=skipped)
    return arr, skipped


def load():
    if not os.path.exists(CACHE):
        return collect()
    z = np.load(CACHE)
    return z["data"], int(z["skipped"])


COL = dict(ratio="#1f6fb4", zone="#2e8b57", cen="#c1440e", ceil="#8a8a8a")

PANEL_DIR = "fig/panels"


def _save(fig, stem, formats=("pdf", "svg", "png"), **kw):
    """Every output is written on a transparent background."""
    for ext in formats:
        fig.savefig(f"{stem}.{ext}", transparent=True, **kw)


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


def panel_cost(ax, x, z2q, z2lo, z2hi, c2q, title=False, standalone=False):
    """(a) two-qubit gates of one oracle pass: centralized against per zone."""
    ax.set_yscale("log")
    ax.axhline(BUDGET_DEFAULT, color=COL["zone"], ls="--", lw=1.0,
               dashes=(4, 3), alpha=0.9)
    ax.text(x[-1], BUDGET_DEFAULT * 1.18, f"partition budget {BUDGET_DEFAULT}",
            ha="right", va="bottom", fontsize=7, color=COL["zone"])

    ax.plot(x, c2q, "s--", color=COL["cen"], ms=5, lw=1.8,
            label="centralized: one circuit for the whole region")
    ax.fill_between(x, z2lo, z2hi, color=_tint(COL["zone"], 0.20), lw=0)
    ax.plot(x, z2q, "o-", color=COL["zone"], ms=5, lw=1.8,
            label="distributed: largest zone circuit")
    ax.annotate(f"×{c2q[-1]/c2q[0]:.0f} across this range", (x[-1], c2q[-1]),
                textcoords="offset points", xytext=(-6, 8), ha="right",
                fontsize=7.5, color=COL["cen"])
    ax.annotate(f"×{z2q[-1]/z2q[0]:.2f}", (x[-1], z2q[-1]),
                textcoords="offset points", xytext=(-6, -13), ha="right",
                fontsize=7.5, color=COL["zone"])
    ax.set_ylabel("two-qubit gates\n(one oracle pass)", fontsize=9.5)
    ax.set_ylim(min(z2q) / 3, max(c2q) * 6)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.25, which="both", lw=0.5)
    if title:
        ax.set_title("Growing the region adds zones, not zone size", fontsize=11)
    if standalone:
        ax.set_xlabel("zones after partitioning", fontsize=10)


def panel_comm(ax, x, tot_kb, per_kb, standalone=False):
    """(b) classical coordination traffic, in total and per zone."""
    ax.set_yscale("log")
    ax.plot(x, tot_kb, "s--", color=COL["cen"], ms=5, lw=1.8,
            label="whole region, one round")
    ax.plot(x, per_kb, "o-", color=COL["zone"], ms=5, lw=1.8,
            label="per zone")
    ax.annotate(f"×{tot_kb[-1]/tot_kb[0]:.1f} for ×{x[-1]/x[0]:.1f} zones",
                (x[-1], tot_kb[-1]), textcoords="offset points",
                xytext=(-6, 8), ha="right", fontsize=7.5, color=COL["cen"])
    ax.annotate(f"×{per_kb[-1]/per_kb[0]:.2f}", (x[-1], per_kb[-1]),
                textcoords="offset points", xytext=(-6, -13), ha="right",
                fontsize=7.5, color=COL["zone"])
    ax.set_ylabel("classical report\n[kB]", fontsize=9.5)
    ax.set_ylim(min(per_kb) / 3, max(tot_kb) * 6)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.25, which="both", lw=0.5)
    if standalone:
        ax.set_xlabel("zones after partitioning", fontsize=10)


def panel_quality(ax, x, ratio, rlo, rhi, n_ues, standalone=False):
    """(c) distributed utility as a fraction of the centralized strict optimum."""
    ax.axhline(100, color="#999999", lw=0.9, ls="-")
    ax.fill_between(x, rlo, rhi, color=_tint(COL["ratio"], 0.22), lw=0)
    ax.plot(x, ratio, "^-", color=COL["ratio"], ms=6, lw=2.0,
            label="distributed utility / centralized strict optimum")
    ax.set_ylim(95, 101)
    ax.set_ylabel("utility vs.\ncentralized optimum  [%]", fontsize=9.5)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{xi:.0f}\n{nu:.0f}" for xi, nu in zip(x, n_ues)],
                       fontsize=8.5)
    ax.set_xlabel("zones after partitioning  /  UEs in the region", fontsize=10)


def plot(arr, skipped):
    cols = dict(g=0, seed=1, n_ue=2, zones=3, ratio=4, zone2q=5, cen2q=6,
                zoneq=7, cenq=8, bdens=9, bits=10, exc=11, feas=12,
                shots_max=13, attempts=14, backed_off=15, ess_min=16,
                beta_min=17)
    gs = np.unique(arr[:, cols["g"]])

    def agg(key):
        m, lo, hi, x = [], [], [], []
        for g in gs:
            v = arr[arr[:, cols["g"]] == g, cols[key]]
            m.append(v.mean())
            lo.append(v.mean() - v.std())
            hi.append(v.mean() + v.std())
            x.append(arr[arr[:, cols["g"]] == g, cols["zones"]].mean())
        return np.array(x), np.array(m), np.array(lo), np.array(hi)

    x, ratio, rlo, rhi = agg("ratio")
    _, z2q, z2lo, z2hi = agg("zone2q")
    _, c2q, _, _ = agg("cen2q")

    n_ues = [arr[arr[:, cols["g"]] == g, cols["n_ue"]].mean() for g in gs]

    _, bits, _, _ = agg("bits")
    tot_kb = bits / 8 / 1024
    per_kb = tot_kb / x

    fig, (ax_c, ax_m, ax_q) = plt.subplots(
        3, 1, figsize=(6.9, 8.0), sharex=True,
        gridspec_kw=dict(height_ratios=[1.15, 1.15, 1]))
    panel_cost(ax_c, x, z2q, z2lo, z2hi, c2q, title=True)
    panel_comm(ax_m, x, tot_kb, per_kb)
    panel_quality(ax_q, x, ratio, rlo, rhi, n_ues)

    feas = arr[:, cols["feas"]].mean() * 100
    fig.text(0.013, 0.055,
             f"{len(arr)} instances, {len(SEEDS)} seeds per size; shading is one "
             f"standard deviation. beta={BETA}, K_z={K_ACCEPT}"
             + (f"; {skipped} globally infeasible instances excluded." if skipped
                else "."),
             fontsize=7, color="#555555")
    fig.text(0.013, 0.022,
             f"Every accepted assignment satisfies the original constraints; "
             f"{feas:.0f}% of runs closed on a strictly feasible global assignment.",
             fontsize=7, color="#555555")
    fig.subplots_adjust(left=0.135, right=0.98, top=0.945, bottom=0.135,
                        hspace=0.14)
    os.makedirs("fig", exist_ok=True)
    _save(fig, "fig/fig_scaling", formats=("pdf", "png"), dpi=200)
    plt.close(fig)
    print("wrote fig/fig_scaling.pdf and .png")

    # ---- the same panels again, standalone and separately editable -------
    os.makedirs(PANEL_DIR, exist_ok=True)
    # clear this script's own panels first: names carry zone and UE indices,
    # so a changed instance would otherwise leave stale files behind
    for old in glob.glob(f"{PANEL_DIR}/scaling_*"):
        os.remove(old)
    specs = [
        ("scaling_a_circuit_cost",
         lambda ax: panel_cost(ax, x, z2q, z2lo, z2hi, c2q, standalone=True)),
        ("scaling_b_coordination_cost",
         lambda ax: panel_comm(ax, x, tot_kb, per_kb, standalone=True)),
        ("scaling_c_utility_ratio",
         lambda ax: panel_quality(ax, x, ratio, rlo, rhi, n_ues, standalone=True)),
    ]
    for stem, draw in specs:
        f, a = plt.subplots(figsize=(5.4, 3.4))
        draw(a)
        f.tight_layout()
        _save(f, f"{PANEL_DIR}/{stem}", dpi=200)
        plt.close(f)
    print(f"wrote {len(specs)} standalone panels to {PANEL_DIR}/ (pdf + svg + png)")

    # ---- console summary -------------------------------------------------
    print("\n zones   UEs   ratio%      max zone 2q     centralized 2q   |B|%")
    for i, g in enumerate(gs):
        m = arr[:, cols["g"]] == g
        print(f" {x[i]:5.1f} {arr[m, cols['n_ue']].mean():5.0f}  "
              f"{ratio[i]:6.2f}+-{arr[m, cols['ratio']].std():4.2f}   "
              f"{z2q[i]:8.0f}         {c2q[i]:10.0f}   "
              f"{arr[m, cols['bdens']].mean():5.1f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--recollect", action="store_true")
    a = ap.parse_args()
    if a.recollect and os.path.exists(CACHE):
        os.remove(CACHE)
    arr, skipped = load()
    plot(arr, skipped)
