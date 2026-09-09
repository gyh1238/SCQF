"""
What the execution budget costs, and what it buys (Sec. V-E, last paragraph).
============================================================================
A zone that runs at the target exponent beta pays K_z / P_z(k) executions for
its report, and in the tightest zones that number is not payable.  Sec. IV-C's
remedy is to execute at the largest beta_z <= beta the budget allows and to
reconstruct the target exponent afterwards by reweighting, Eq. (A.10).

The trade is then explicit and worth measuring rather than asserting:

    uncapped   every zone at beta, whatever it costs
    capped     every zone at 10^4 executions, backing off beta_z as needed

Reported for both: the executions the most expensive single zone demands, and
the utility the protocol reaches.  Reconstruction is unbiased, so the cap
should cost effective sample size and not accuracy; that is the claim this
module checks.

Run:  python backoff_study.py
"""

import numpy as np

from model_instance import make_instance, utility_scale
from model_partition import partition_aps
from proto_coordination import run_protocol
from model_reference import solve_centralized
from proto_zone import (K_ROUNDS, amplification_rounds, build_zone,
                       choose_execution_exponent)

BETA = 1.5
K_ACCEPT = 200
CAP = 10_000
UNCAPPED = np.inf
G_VALUES = (4, 5, 6, 7)
SEEDS = range(4)


def instances():
    for g in G_VALUES:
        for s in SEEDS:
            inst = make_instance(g=g, seed=s)
            opt, _, ok = solve_centralized(inst)
            if ok:
                yield g, s, inst, partition_aps(inst), opt


def demand_at_target(inst, part, ubar, k=K_ROUNDS):
    """Executions each zone would need at the target exponent, uncapped."""
    out = []
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue == 0:
            continue
        _, mu, shots, _ = choose_execution_exponent(
            zone, BETA, ubar, K_ACCEPT, UNCAPPED, k=k)
        out.append((shots, mu, zone.n_ue))
    return out


def main():
    print("=" * 74)
    print(f"  Execution budget vs utility, k = {K_ROUNDS}, K_z = {K_ACCEPT}")
    print("=" * 74)

    peak_unc, peak_cap, r_unc, r_cap, rounds = [], [], [], [], []
    for g, s, inst, part, opt in instances():
        ubar = utility_scale(inst)
        d = demand_at_target(inst, part, ubar)
        peak_unc.append(max(x[0] for x in d))
        rounds += [amplification_rounds(mu) for _, mu, _ in d]

        for budget, bucket, ratios in ((UNCAPPED, peak_unc, r_unc),
                                       (CAP, peak_cap, r_cap)):
            res = run_protocol(inst, part, beta=BETA, ubar=ubar,
                               k_accept=K_ACCEPT, shot_budget=budget,
                               rng=np.random.default_rng(1000 + s))
            ratios.append(100.0 * res["utility"] / opt)
            if budget is CAP:
                bucket.append(res["shots_max"])

    peak_unc = np.array(peak_unc)
    print(f"\n  instances                {len(r_unc)}")
    print(f"  zones                    {len(rounds)}")
    print(f"\n  executions at the target exponent, per zone")
    print(f"    range over zones       {min(peak_unc):,.0f} to {max(peak_unc):,.0f}")
    print(f"    peak single zone       {max(peak_unc):,.3g}")
    print(f"    capped at {CAP:,}       reduction x{max(peak_unc) / CAP:,.0f}")
    print(f"\n  rounds maximizing P_z(k), over all zones")
    print(f"    range                  {min(rounds)} to {max(rounds)}")
    print(f"\n  utility, mean over instances")
    print(f"    uncapped               {np.mean(r_unc):.2f}% of J*")
    print(f"    capped at {CAP:,}       {np.mean(r_cap):.2f}% of J*")
    d = np.array(r_unc) - np.array(r_cap)
    print(f"    loss                   {d.mean():.2f} pp on average, "
          f"{d.max():.2f} pp at worst")


if __name__ == "__main__":
    main()
