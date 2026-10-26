"""Plot the zone-local classical--quantum crossover for joint AP--RB assignment.

The plotted instance family uses two joint candidates j=(a, r) per UE, so
D_z=2N_z.  Proposed, GAS, and QTG use the same capacity-aware three-bit
resource-counter profile.  The enhanced Proposed schedule additionally uses
four resource counters, UE-parallel match/utility layers, and a 0.75 arithmetic
depth factor.  The same optimized arithmetic primitive is applied to GAS and
QTG, while their algorithm-specific sequential structures remain.

Quantum curves are modeled one-pass logical latencies.  The classical curve is
the representative scenario 0.400 N_z^3 ns, a rounded interior coefficient
within the previously considered 0.058--0.740 interval.
"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "unified-assignment-matplotlib-cache"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
import numpy as np


TAU_NS = 12.5
CLASSICAL_COEFFICIENT_NS = 0.400
CANDIDATES_PER_UE = 2.0
RESOURCE_INCIDENCE = 2.0
CAPACITY_CONDITIONS_PER_UE = 1.0
PARALLEL_COUNTERS = 4.0
RESOURCE_COUNTER_BITS = 3.0
ARITHMETIC_DEPTH_FACTOR = 0.75


def profile(zone_ues: np.ndarray) -> tuple[np.ndarray, ...]:
    """Continuous bounded-density profile used by all quantum curves."""
    n = np.asarray(zone_ues, dtype=float)
    d = CANDIDATES_PER_UE * n
    c = CAPACITY_CONDITIONS_PER_UE * n
    ell = np.full_like(n, np.log2(CANDIDATES_PER_UE))
    b = np.full_like(n, RESOURCE_COUNTER_BITS)
    # GAS and QTG still accumulate a zone-wide numeric objective/profit.
    p = np.log2(n + 1.0)
    return n, d, c, ell, b, p


def flag_depth(n: np.ndarray, c: np.ndarray) -> np.ndarray:
    return 2.0 * np.log2(n + c + 1.0)


def proposed_depth(
    zone_ues: np.ndarray,
    counters: float = 1.0,
    ue_parallel: bool = False,
) -> np.ndarray:
    """One Proposed pass; only accumulation is divided by the counter count."""
    n, d, c, ell, b, _ = profile(zone_ues)
    if ue_parallel:
        matching = np.full_like(n, 2.0 * CANDIDATES_PER_UE * ell[0])
        utility_rotations = np.full_like(n, 2.0 * CANDIDATES_PER_UE)
    else:
        matching = 2.0 * d * ell
        utility_rotations = 2.0 * d
    accumulation = 2.0 * RESOURCE_INCIDENCE * d * b / counters
    comparisons = 2.0 * c * b
    arithmetic = ARITHMETIC_DEPTH_FACTOR * (accumulation + comparisons)
    return matching + utility_rotations + arithmetic + flag_depth(n, c)


def gas_depth(zone_ues: np.ndarray) -> np.ndarray:
    """One GAS threshold-oracle pass for the joint candidate graph."""
    n, d, c, ell, b, p = profile(zone_ues)
    matching = 2.0 * d * ell
    constraint_accumulation = 2.0 * RESOURCE_INCIDENCE * d * b
    objective_accumulation = 2.0 * d * p
    comparisons = 2.0 * (c + 1.0) * b
    arithmetic = ARITHMETIC_DEPTH_FACTOR * (
        constraint_accumulation + objective_accumulation + comparisons
    )
    return matching + arithmetic + flag_depth(n, c)


def qtg_depth(zone_ues: np.ndarray) -> np.ndarray:
    """One search pass of the joint-capacity QTG extension."""
    n, d, c, ell, b, p = profile(zone_ues)
    categorical_checks = 2.0 * d * ell
    capacity_checks_updates_and_profit = (
        2.0 * d * (2.0 * RESOURCE_INCIDENCE * b + p)
    )
    arithmetic = ARITHMETIC_DEPTH_FACTOR * capacity_checks_updates_and_profit
    return categorical_checks + arithmetic + flag_depth(n, c)


def quantum_seconds(depth: np.ndarray) -> np.ndarray:
    return depth * TAU_NS * 1.0e-9


def classical_seconds(zone_ues: np.ndarray) -> np.ndarray:
    n = np.asarray(zone_ues, dtype=float)
    return CLASSICAL_COEFFICIENT_NS * n**3 * 1.0e-9


def first_crossover(
    n_values: np.ndarray,
    quantum: np.ndarray,
    classical: np.ndarray,
) -> float:
    """Interpolate the first quantum/classical equality on log-log axes."""
    indices = np.flatnonzero(quantum <= classical)
    if not indices.size:
        raise RuntimeError("no crossover in plotted domain")
    index = int(indices[0])
    if index == 0:
        return float(n_values[0])

    x0, x1 = np.log(n_values[index - 1 : index + 1])
    ratio = np.log(
        quantum[index - 1 : index + 1]
        / classical[index - 1 : index + 1]
    )
    return float(np.exp(x0 - ratio[0] * (x1 - x0) / (ratio[1] - ratio[0])))


def draw(output_dir: Path) -> tuple[Path, Path]:
    # Continuous relaxed widths produce smooth trend curves rather than
    # artificial power-of-two stair steps.
    n_values = np.geomspace(3.0, 1000.0, 3000)

    proposed_serial = quantum_seconds(
        proposed_depth(n_values, counters=1.0, ue_parallel=False)
    )
    proposed_enhanced = quantum_seconds(
        proposed_depth(
            n_values,
            counters=PARALLEL_COUNTERS,
            ue_parallel=True,
        )
    )
    gas = quantum_seconds(gas_depth(n_values))
    qtg = quantum_seconds(qtg_depth(n_values))
    classical = classical_seconds(n_values)

    crossovers = {
        "Proposed, enhanced": first_crossover(
            n_values, proposed_enhanced, classical
        ),
        "Proposed, serial": first_crossover(n_values, proposed_serial, classical),
        "GAS": first_crossover(n_values, gas, classical),
        "QTG": first_crossover(n_values, qtg, classical),
    }

    # Drawn at the one-column width it is printed at (IEEE, 3.5 in), so every
    # size here is the size on paper.  The axes are cropped to the decade in
    # which the four crossovers fall, which is what the figure is read for.
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 0.8,
        }
    )

    fig, ax = plt.subplots(figsize=(3.5, 2.75))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(8.0, 80.0)
    ax.set_ylim(2.0e-7, 2.5e-3)      # headroom: the legend sits above every curve

    colors = {
        "Proposed, serial": "#1769AA",
        "Proposed, enhanced": "#2A9D8F",
        "GAS": "#E67E22",
        "QTG": "#7B4AB5",
        "Classical": "#D62728",
    }

    ax.plot(n_values, proposed_serial, color=colors["Proposed, serial"],
            linewidth=1.5, label="Proposed serial", zorder=4)
    ax.plot(n_values, proposed_enhanced, color=colors["Proposed, enhanced"],
            linewidth=1.5, linestyle=(0, (5, 2)),
            label=r"Proposed enhanced ($Q^{\mathrm{cnt}}=3$, $R_{\mathrm{cnt}}=4$)",
            zorder=4)
    ax.plot(n_values, gas, color=colors["GAS"], linewidth=1.5, label="GAS",
            zorder=3)
    ax.plot(n_values, qtg, color=colors["QTG"], linewidth=1.5, label="QTG",
            zorder=3)
    ax.plot(n_values, classical, color=colors["Classical"], linewidth=1.8,
            label=r"Classical: $0.400N_z^3$ ns", zorder=2)

    # Each crossover is marked on its curve and its N_z written upright on its
    # dashed line, in the empty lower part of the axes.  The three right-hand
    # lines are closer together than a label is wide, so their labels step up
    # one above the other instead of sitting side by side.
    label_y = {
        "Proposed, enhanced": 3.6e-7,
        "Proposed, serial": 3.6e-7,
        "GAS": 7.5e-7,
        "QTG": 1.55e-6,
    }
    for name, crossing in crossovers.items():
        color = colors[name]
        y_cross = float(classical_seconds(np.array([crossing]))[0])
        ax.axvline(crossing, color=color, linewidth=0.8,
                   linestyle=(0, (3, 3)), alpha=0.8, zorder=1)
        ax.scatter([crossing], [y_cross], s=22, facecolor="white",
                   edgecolor=color, linewidth=1.2, zorder=7)
        ax.text(crossing, label_y[name], f"{crossing:.1f}", ha="center",
                va="center", fontsize=8, color=color, fontweight="bold",
                zorder=8, bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                    ec="none"))

    ax.set_xlabel(r"Zone size $N_z$ (UEs)")
    ax.set_ylabel("Modeled logical latency (s)")
    ticks = [8, 10, 20, 30, 40, 50, 80]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

    top_axis = ax.secondary_xaxis(
        "top",
        functions=(
            lambda n: CANDIDATES_PER_UE * n,
            lambda d: d / CANDIDATES_PER_UE,
        ),
    )
    top_axis.set_xlabel(r"Candidate-edge count $D_z=2N_z$", labelpad=3)
    top_axis.xaxis.set_major_locator(FixedLocator([20, 40, 60, 100, 160]))
    top_axis.xaxis.set_major_formatter(ScalarFormatter())
    top_axis.xaxis.set_minor_formatter(NullFormatter())
    top_axis.tick_params(labelsize=7.5)

    ax.grid(which="major", color="#A9A9A9", linestyle=":", linewidth=0.6,
            alpha=0.8)
    ax.grid(which="minor", axis="y", color="#D5D5D5", linestyle=":",
            linewidth=0.4, alpha=0.6)

    legend = ax.legend(loc="upper left", frameon=True, framealpha=0.95,
                       borderpad=0.35, handlelength=2.0, labelspacing=0.25,
                       borderaxespad=0.3)
    legend.get_frame().set_edgecolor("#999999")
    legend.get_frame().set_linewidth(0.6)

    fig.tight_layout(pad=0.2)

    png_path = output_dir / "qtg_assignment_runtime_comparison.png"
    pdf_path = output_dir / "qtg_assignment_runtime_comparison.pdf"
    # transparent like the other figures, so all of them sit on the page alike
    fig.savefig(png_path, dpi=400, bbox_inches="tight", pad_inches=0.02,
                transparent=True)
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.02,
                transparent=True)
    plt.close(fig)

    for name, crossing in crossovers.items():
        print(f"{name}: N_z={crossing:.3f}, D_z={2.0 * crossing:.3f}")
    print(f"saved: {png_path}")
    print(f"saved: {pdf_path}")
    return png_path, pdf_path


if __name__ == "__main__":
    draw(Path(__file__).resolve().parent)
