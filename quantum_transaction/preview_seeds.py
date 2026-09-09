"""
Contact sheet of candidate snapshot regions.
============================================
The snapshot figure shows one instance, and which one is a presentation
choice as much as a statistical one: every seed is a valid draw from the
same generator, but they do not all read equally well as a picture.

This renders the partition panel for a range of seeds side by side so the
region can be chosen by eye, and prints the numbers that matter for the
choice: how many APs the ground actually accepted, how evenly the
partitioner split the region, whether it left many single-AP zones, and how
heavy the boundary is.

Instances come from `fig_snapshot.build_instance`, so the preview is
drawn on whatever geometry the figure uses -- with `snap.GEO` set, the same
window of campus, over the same basemap.

Set the chosen value as `SEED` in `fig_snapshot.py`.

Usage:  python preview_seeds.py [n_seeds]
"""

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model_instance import utility_scale
from model_partition import partition_aps
from proto_coordination import run_protocol
from model_reference import solve_centralized
import fig_snapshot as snap


def evaluate(seed):
    inst = snap.build_instance(seed)
    opt, _, ok = solve_centralized(inst)
    if not ok:
        return None
    part = partition_aps(inst)
    res = run_protocol(inst, part, beta=snap.BETA, ubar=utility_scale(inst),
                       k_accept=snap.K_ACCEPT, shot_budget=snap.SHOT_BUDGET,
                       rng=np.random.default_rng(1000 + seed), collect=True)
    active = [r for r in res["reports"] if r["zone"].n_ue > 0]
    n_ap = [len(z) for z in part.zones if z]
    n_ue = [r["zone"].n_ue for r in active]
    return dict(seed=seed, inst=inst, part=part, res=res, active=active,
                ratio=100.0 * res["utility"] / opt, n_ap=inst.n_ap,
                zones=len(active), singles=sum(1 for a in n_ap if a == 1),
                ap_spread=float(np.std(n_ap)), ue_spread=float(np.std(n_ue)),
                bdens=100.0 * part.boundary_density)


def main(n_seeds=8):
    scored = [(s, evaluate(s)) for s in range(n_seeds)]
    cands = [c for _, c in scored if c]
    dropped = [s for s, c in scored if c is None]
    if dropped:
        # On the campus a draw can be genuinely infeasible -- UEs bunch onto
        # open ground faster than the APs facing it can admit them -- and the
        # centralized reference then has no optimum to divide by.  Say which
        # seeds went, so a thin contact sheet does not look like a bug.
        print(f"no centralized optimum for seeds {dropped}; "
              f"{len(cands)}/{n_seeds} shown\n")
    print(f"{'seed':>4} {'APs':>4} {'zones':>6} {'1-AP':>5} {'sd(APs)':>8} "
          f"{'sd(N_z)':>8} {'|B|%':>6} {'ratio%':>7}")
    for c in cands:
        print(f"{c['seed']:>4} {c['n_ap']:>4} {c['zones']:>6} {c['singles']:>5} "
              f"{c['ap_spread']:>8.2f} {c['ue_spread']:>8.2f} "
              f"{c['bdens']:>6.0f} {c['ratio']:>7.2f}")

    ncol = 4
    nrow = int(np.ceil(len(cands) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 4.3 * nrow))
    cmap = plt.get_cmap("tab20")
    for ax, c in zip(np.ravel(axes), cands):
        zcol = {z: cmap(i % 20) for i, z in
                enumerate(sorted(r["zone"].idx for r in c["active"]))}
        snap.panel_map(ax, c["inst"], c["part"], c["active"], zcol)
        ax.get_legend().remove()
        ax.set_title(f"seed {c['seed']}: {c['zones']} zones, "
                     f"{c['singles']} single-AP, |B|={c['bdens']:.0f}%",
                     fontsize=10, pad=5)
    for ax in np.ravel(axes)[len(cands):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig("fig/seed_preview.png", dpi=160, transparent=True)
    print("\nwrote fig/seed_preview.png")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
