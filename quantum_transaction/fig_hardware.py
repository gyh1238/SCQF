"""
How much of the utility weighting survives on hardware, against circuit size.
============================================================================
Every inter-cell run in `hw_runs/` that returned counts, scored by the same
code as the tables (`noise_oracle.law_from_counts`).  Each run gets a grey bar
for what its own circuit could show: one end is the noise-free value, the other
the value of an output carrying no information, which after the acceptance test
is uniform over the feasible set.  The marker says where on that range the run
landed, so runs on circuits with different feasible sets can share an axis.

Run:
    python fig_hardware.py        # writes fig/fig_hardware.{pdf,png}
"""

import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
import numpy as np

import circuit_inter as inter
from certify_rotation import gibbs
from noise_oracle import INTER_UTILITY, law_from_counts

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_W = 3.5           # printed at one column, so every size is the size on paper
FS_LABEL = 8.0
FS_TICK = 7.5
FS_NOTE = 7.0
FS_LETTER = 8.5

COLORS = {"ibm_boston": "#1769AA", "ibm_kingston": "#E67E22", "IonQ Forte": "#7B4AB5"}
BAR = "#C9CDD3"
N_BOOT = 2000


def accepted(counts, n, feas):
    """Counts of the accepted assignments, by decoded code."""
    out = {}
    for bits, c in counts.items():
        bits = bits.replace(" ", "")
        if "1" in bits[::-1][n:2 * n]:
            continue
        code = bits[-n:]
        if code in feas:
            out[code] = out.get(code, 0) + c
    return out


def load():
    runs = []
    rng = np.random.default_rng(0)
    for path in sorted(glob.glob(os.path.join(HERE, "hw_runs", "*.json"))):
        with open(path, encoding="utf-8") as f:
            r = json.load(f)
        if not r.get("counts") or r.get("simulator") or "simulator" in r["backend"]:
            continue
        if not r["case"].startswith("inter"):
            continue
        w, cap = r.get("demands") or [1, 1, 1, 1], r.get("cap") or 2
        u = INTER_UTILITY[:len(w)]
        n = u.shape[0]                  # the table has four rows: at most four UEs
        feas = {"".join(str(b) for b in reversed(bits)): s
                for bits, s in inter.enumerate_feasible(u, w, cap)}
        meta = dict(n_assign=n, n_cost=n, lam=r["lam"])
        m = law_from_counts(r["counts"], meta, feas)
        assert abs(m["tvd"] - r["metrics"]["tvd"]) < 1e-9, path

        ref = gibbs(feas, r["lam"])
        opt = max(ref, key=ref.get)
        acc = accepted(r["counts"], n, feas)
        n_acc = sum(acc.values())
        keys = sorted(ref)
        p_ref = np.array([ref[k] for k in keys])
        p_hat = np.array([acc.get(k, 0) / n_acc for k in keys])
        boot = rng.multinomial(n_acc, p_hat, N_BOOT) / n_acc
        i_opt = keys.index(opt)
        uniform = 0.5 * np.abs(p_ref - 1 / len(feas)).sum()

        ionq = "ionq" in r["backend"] or "forte" in r["backend"]
        g = r.get("two_qubit") if ionq else r.get("routed_2q")
        dd = r.get("options", {}).get("dynamical_decoupling", True)
        lo_u = np.exp(-r["lam"] * (u.max(axis=1)[:, None] - u))
        mu = sum(np.prod([lo_u[i, int(c[::-1][i])] for i in range(n)])
                 for c in feas) / 2 ** n
        runs.append(dict(
            device="IonQ Forte" if ionq else r["backend"], k=r["k"], gates=g, dd=dd,
            n=n, F=len(feas), n_acc=n_acc, demands=w, seed=r.get("options", {}).get("seed_transpiler"),
            accept=m["accept"], accept_random=2.0 ** (-2 * n) * len(feas),
            mu=mu, p1=np.sin(3 * np.arcsin(np.sqrt(mu))) ** 2,
            tvd=m["tvd"], tvd_err=(0.5 * np.abs(boot - p_ref).sum(axis=1)).std(),
            tvd_uniform=uniform,
            popt=p_hat[i_opt], popt_err=boot[:, i_opt].std(),
            popt_ideal=p_ref[i_opt], popt_uniform=1 / len(feas)))
    return runs


def spread(runs, ratio=1.06):
    """Runs of nearly equal size are set side by side around their true count."""
    runs.sort(key=lambda r: (r["gates"], r["device"], tuple(r["demands"])))
    clusters = [[runs[0]]]
    for r in runs[1:]:
        if r["gates"] / clusters[-1][-1]["gates"] < ratio:
            clusters[-1].append(r)
        else:
            clusters.append([r])
    for rs in clusters:
        mid = np.exp(np.mean([np.log(r["gates"]) for r in rs]))
        for j, r in enumerate(rs):
            r["x"] = mid * ratio ** (j - (len(rs) - 1) / 2)


def panel(ax, runs, val, err, lo, hi):
    for r in runs:
        ax.plot([r["x"], r["x"]], [r[lo], r[hi]], color=BAR, lw=2.4,
                solid_capstyle="butt", zorder=1)
    for r in runs:
        c = COLORS[r["device"]]
        ax.errorbar(r["x"], r[val], yerr=r[err], color=c, lw=0.7, zorder=3,
                    marker="o", ms=3.6, mfc=c if r["dd"] else "white", mec=c,
                    mew=0.8, capsize=0)
    ax.tick_params(labelsize=FS_TICK, length=2.5, width=0.6, pad=2)
    ax.grid(axis="y", color="#D9D9D9", lw=0.5, ls=":")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    runs = [r for r in load() if r["gates"] is not None]
    for r in runs:
        r["zero"] = 0.0
    k0 = sorted((r for r in runs if r["k"] == 0),
                key=lambda r: (r["gates"], r["device"], tuple(r["demands"])))
    k1 = [r for r in runs if r["k"] == 1]
    for j, r in enumerate(k0):
        r["x"] = j
    spread(k1)
    for r in sorted(runs, key=lambda r: (r["device"], r["k"], r["gates"])):
        print(f"{r['device']:<12} k={r['k']} gates {r['gates']:>5} dd={r['dd']!s:<5} "
              f"seed={r['seed']} w={r['demands']} |F|={r['F']:>2} n_acc={r['n_acc']:>5}  "
              f"acc {r['accept']:.3f} (rand {r['accept_random']:.3f}, mu {r['mu']:.3f}, "
              f"P1 {r['p1']:.3f})  TVD {r['tvd']:.3f}+-{r['tvd_err']:.3f} "
              f"(unif {r['tvd_uniform']:.3f})  Popt {r['popt']:.3f}+-{r['popt_err']:.3f} "
              f"(ideal {r['popt_ideal']:.3f}, unif {r['popt_uniform']:.3f})")

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                         "mathtext.fontset": "dejavusans"})
    fig, ax = plt.subplots(2, 2, figsize=(FIG_W, 3.3), sharey="row", sharex="col",
                           gridspec_kw=dict(width_ratios=[1, 2.5], wspace=0.06))
    rows = [("tvd", "tvd_err", "zero", "tvd_uniform", (0, 0.56), r"TVD from $p_z$"),
            ("popt", "popt_err", "popt_uniform", "popt_ideal", (0, 0.8),
             r"$P(\mathrm{opt})$")]
    for (val, err, lo, hi, ylim, label), (left, right) in zip(rows, ax):
        panel(left, k0, val, err, lo, hi)
        panel(right, k1, val, err, lo, hi)
        left.set_ylim(*ylim)
        left.set_ylabel(label, fontsize=FS_LABEL, labelpad=2)
        left.set_xlim(-0.8, len(k0) - 0.2)
        left.set_xticks([])
        right.set_xscale("log")
        right.set_xlim(160.0, 2700.0)
        right.xaxis.set_major_locator(FixedLocator([200, 500, 1000, 2000]))
        right.xaxis.set_major_formatter(ScalarFormatter())
        right.xaxis.set_minor_formatter(NullFormatter())
        right.tick_params(axis="y", length=0)
        right.spines["left"].set_visible(False)
    ax[1, 0].set_xlabel("$k=0$\n4 to 8 gates", fontsize=FS_LABEL, labelpad=2)
    ax[1, 1].set_xlabel("$k=1$\ntwo-qubit gates", fontsize=FS_LABEL, labelpad=2)
    fig.align_xlabels(ax[1, :])
    ax[0, 0].text(0.04, 0.97, "(a)", transform=ax[0, 0].transAxes,
                  fontsize=FS_LETTER, va="top")
    ax[1, 0].text(0.04, 0.03, "(b)", transform=ax[1, 0].transAxes,
                  fontsize=FS_LETTER, va="bottom")

    nodd = next(r for r in runs if not r["dd"])
    ax[0, 1].annotate("no DD", (nodd["x"], nodd["tvd"]), xytext=(5, 0), fontsize=FS_NOTE,
                      textcoords="offset points", ha="left", va="center")
    ax[1, 1].annotate("no DD", (nodd["x"], nodd["popt"]), xytext=(5, 0), fontsize=FS_NOTE,
                      textcoords="offset points", ha="left", va="center")

    handles = [Line2D([], [], color=c, marker="o", ls="", ms=3.6, label=d)
               for d, c in COLORS.items()]
    handles.append(Line2D([], [], color=BAR, lw=2.4, label="noise-free to random"))
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=6.6,
               frameon=False, handletextpad=0.35, columnspacing=0.7, handlelength=1.0,
               borderaxespad=0.0, bbox_to_anchor=(0.54, 1.0))
    fig.subplots_adjust(left=0.125, right=0.985, top=0.925, bottom=0.135,
                        hspace=0.12, wspace=0.06)
    os.makedirs(os.path.join(HERE, "fig"), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig", f"fig_hardware.{ext}"), dpi=400,
                    bbox_inches="tight", pad_inches=0.02, transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
