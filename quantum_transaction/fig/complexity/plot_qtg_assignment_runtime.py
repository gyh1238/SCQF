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

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.4,
            "axes.labelsize": 12,
            "axes.titlesize": 13.2,
            "legend.fontsize": 9.2,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.0,
        }
    )

    fig, ax = plt.subplots(figsize=(10.0, 6.1))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3.0, 1000.0)
    ax.set_ylim(5.0e-9, 6.0e-1)

    colors = {
        "Proposed, serial": "#1769AA",
        "Proposed, enhanced": "#2A9D8F",
        "GAS": "#E67E22",
        "QTG": "#7B4AB5",
        "Classical": "#D62728",
    }

    ax.plot(
        n_values,
        proposed_serial,
        color=colors["Proposed, serial"],
        linewidth=2.7,
        label="Proposed unified, serial accumulation",
        zorder=4,
    )
    ax.plot(
        n_values,
        proposed_enhanced,
        color=colors["Proposed, enhanced"],
        linewidth=2.5,
        linestyle=(0, (7, 3)),
        label=r"Proposed enhanced: capacity-aware, $s=4$",
        zorder=4,
    )
    ax.plot(
        n_values,
        gas,
        color=colors["GAS"],
        linewidth=2.5,
        label="GAS, one threshold-oracle pass",
        zorder=3,
    )
    ax.plot(
        n_values,
        qtg,
        color=colors["QTG"],
        linewidth=2.5,
        label="QTG, one joint-capacity search pass",
        zorder=3,
    )
    ax.plot(
        n_values,
        classical,
        color=colors["Classical"],
        linewidth=2.9,
        label=r"Classical scenario: $0.400N_z^3$ ns",
        zorder=2,
    )

    enhanced_cross = crossovers["Proposed, enhanced"]
    favorable = n_values >= enhanced_cross
    ax.fill_between(
        n_values[favorable],
        proposed_enhanced[favorable],
        classical[favorable],
        color=colors["Proposed, enhanced"],
        alpha=0.065,
        linewidth=0,
        zorder=0,
    )

    label_offsets = {
        "Proposed, enhanced": (0.82, 0.49),
        "Proposed, serial": (0.91, 2.30),
        "GAS": (1.03, 0.42),
        "QTG": (1.10, 2.35),
    }

    for name, crossing in crossovers.items():
        color = colors[name]
        y_cross = float(classical_seconds(np.array([crossing]))[0])
        ax.axvline(
            crossing,
            color=color,
            linewidth=1.1,
            linestyle=(0, (4, 4)),
            alpha=0.72,
            zorder=1,
        )
        ax.scatter(
            [crossing],
            [y_cross],
            s=46,
            facecolor="white",
            edgecolor=color,
            linewidth=1.8,
            zorder=7,
        )
        x_factor, y_factor = label_offsets[name]
        ax.text(
            crossing * x_factor,
            y_cross * y_factor,
            rf"$N_z\!\approx\!{crossing:.0f}$",
            color=color,
            fontsize=9.0,
            fontweight="bold",
            ha="center",
            va="center",
        )

    ax.annotate(
        "Capacity-aware enhanced schedule",
        xy=(enhanced_cross, float(classical_seconds(np.array([enhanced_cross]))[0])),
        xytext=(7.5, 2.5e-4),
        color=colors["Proposed, enhanced"],
        fontsize=9.8,
        fontweight="bold",
        arrowprops=dict(
            arrowstyle="->",
            color=colors["Proposed, enhanced"],
            linewidth=1.5,
        ),
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor=colors["Proposed, enhanced"],
            alpha=0.94,
        ),
    )

    crossover_text = (
        "Modeled crossover $N_z$\n"
        f"Proposed enhanced: {crossovers['Proposed, enhanced']:.1f}\n"
        f"Proposed serial: {crossovers['Proposed, serial']:.1f}\n"
        f"GAS: {crossovers['GAS']:.1f}\n"
        f"QTG: {crossovers['QTG']:.1f}"
    )
    ax.text(
        0.975,
        0.76,
        crossover_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.2,
        color="#3F3F3F",
        linespacing=1.35,
        bbox=dict(
            boxstyle="round,pad=0.42",
            facecolor="white",
            edgecolor="#999999",
            alpha=0.94,
        ),
    )

    ax.text(
        0.022,
        0.965,
        r"Capacity-aware profile: $b_m=3$, arithmetic depth $\times0.75$; "
        r"enhanced Proposed also uses $s=4$ and UE-parallel layers",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.0,
        color="#4D4D4D",
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="#A0A0A0",
            alpha=0.94,
        ),
    )

    ax.set_title("Zone-local classical–quantum crossover")
    ax.set_xlabel(r"Zone size $N_z$ (UEs)")
    ax.set_ylabel("Modeled runtime (s)")

    top_axis = ax.secondary_xaxis(
        "top",
        functions=(
            lambda n: CANDIDATES_PER_UE * n,
            lambda d: d / CANDIDATES_PER_UE,
        ),
    )
    top_axis.set_xlabel(r"Candidate-edge count $D_z=2N_z$")

    ax.grid(
        which="major",
        color="#A9A9A9",
        linestyle=":",
        linewidth=0.85,
        alpha=0.78,
    )
    ax.grid(
        which="minor",
        color="#D5D5D5",
        linestyle=":",
        linewidth=0.55,
        alpha=0.62,
    )

    legend = ax.legend(
        loc="lower right",
        frameon=True,
        fancybox=True,
        framealpha=0.96,
        borderpad=0.65,
    )
    legend.get_frame().set_edgecolor("#999999")

    fig.text(
        0.5,
        0.026,
        r"Quantum curves: $\tau=12.5$ ns/layer. Classical coefficient "
        r"$0.400$: rounded representative value within the $0.058$--$0.740$ interval.",
        ha="center",
        va="bottom",
        fontsize=9.0,
        color="#4F4F4F",
    )
    fig.text(
        0.5,
        0.009,
        "GAS threshold repetitions, QTG amplification repetitions, accepted draws, shots, and coordination are excluded.",
        ha="center",
        va="bottom",
        fontsize=8.8,
        color="#5A5A5A",
    )
    fig.tight_layout(rect=(0.015, 0.072, 0.995, 0.98))

    png_path = output_dir / "qtg_assignment_runtime_comparison.png"
    pdf_path = output_dir / "qtg_assignment_runtime_comparison.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    for name, crossing in crossovers.items():
        print(f"{name}: N_z={crossing:.3f}, D_z={2.0 * crossing:.3f}")
    print(f"saved: {png_path}")
    print(f"saved: {pdf_path}")
    return png_path, pdf_path


if __name__ == "__main__":
    draw(Path(__file__).resolve().parent)
