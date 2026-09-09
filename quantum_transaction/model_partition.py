"""
Budget-constrained zone partitioning (the missing upstream of Sec. III-C).
==========================================================================
The manuscript starts from "Partition the APs into zones Z" and takes the
partition as given, while Sec. V invokes the resulting bounded-density
regime as evidence.  This module produces that partition.

Because a UE is a boundary UE exactly when its candidate RBs are owned by
APs in different zones, partitioning the APs *is* the choice of boundary
set.  Building the AP co-coverage graph

    w(a, a') = #{ i : V_i touches both a and a' }

turns "minimize the number of boundary UEs" into a graph cut, and the
hardware limit enters as a capacity constraint on each part:

    minimize   |B| = #{ i : i appears in >= 2 zones }
    subject to n2q(zone) <= budget         for every zone

The budget is a two-qubit gate count, not a qubit count, because that is
what `HARDWARE_LIMITS.md` measured to be binding.

Method: greedy co-coverage region growing under the budget, followed by
boundary-AP refinement passes that move single APs between adjacent zones
whenever the cut strictly decreases.
"""

from dataclasses import dataclass

import numpy as np

from model_cost import zone_2q_cost, BUDGET_DEFAULT


@dataclass
class Partition:
    zone_of: np.ndarray        # (n_ap,) zone index of each AP
    zones: list                # zones[z] = sorted list of AP indices
    ue_zones: list             # ue_zones[i] = sorted list of zones concerning UE i
    boundary: np.ndarray       # indices of boundary UEs
    stats: list                # per-zone dict from zone_2q_cost(detail=True)
    budget: int

    @property
    def n_zones(self):
        return len(self.zones)

    @property
    def boundary_density(self):
        return len(self.boundary) / len(self.ue_zones)

    def summary(self):
        n2q = [s["n2q"] for s in self.stats]
        nue = [s["n_ue"] for s in self.stats]
        return dict(n_zones=self.n_zones, max_2q=max(n2q), med_2q=float(np.median(n2q)),
                    max_ue=max(nue), med_ue=float(np.median(nue)),
                    max_qubits=max(s["n_qubits"] for s in self.stats),
                    n_boundary=len(self.boundary),
                    boundary_density=self.boundary_density)


def ue_owner_aps(inst):
    """ue_aps[i] = unique AP indices owning the candidate RBs of UE i."""
    return [np.unique(inst.rb_owner[c]) for c in inst.cand]


def cocoverage(inst, ue_aps=None):
    """AP co-coverage matrix w(a,a')."""
    if ue_aps is None:
        ue_aps = ue_owner_aps(inst)
    w = np.zeros((inst.n_ap, inst.n_ap))
    for aps in ue_aps:
        for x in range(len(aps)):
            for y in range(x + 1, len(aps)):
                w[aps[x], aps[y]] += 1
                w[aps[y], aps[x]] += 1
    return w


def _cut(zone_of, ue_aps):
    """Number of boundary UEs under the current AP-to-zone map."""
    n = 0
    for aps in ue_aps:
        if len(set(zone_of[a] for a in aps)) > 1:
            n += 1
    return n


def partition_aps(inst, budget=BUDGET_DEFAULT, refine_passes=6, verbose=False):
    """Greedy budget-constrained region growing + boundary refinement."""
    ue_aps = ue_owner_aps(inst)
    w = cocoverage(inst, ue_aps)
    n_ap = inst.n_ap

    zone_of = np.full(n_ap, -1, dtype=int)
    zones = []
    unassigned = set(range(n_ap))

    # ---- greedy growth ---------------------------------------------------
    while unassigned:
        # seed: unassigned AP with the least co-coverage to other unassigned
        # APs, i.e. start from the periphery so that dense cores stay intact
        seed = min(unassigned, key=lambda a: (w[a, list(unassigned)].sum(), a))
        cur = [seed]
        unassigned.discard(seed)
        if zone_2q_cost(inst, cur) > budget:
            # a single AP already exceeds the budget: keep it as its own zone
            zones.append(cur)
            zone_of[seed] = len(zones) - 1
            continue
        while True:
            # candidates ranked by co-coverage with the current zone
            cands = sorted(unassigned, key=lambda a: (-w[a, cur].sum(),
                                                      np.linalg.norm(inst.ap_xy[a] - inst.ap_xy[cur].mean(0))))
            added = False
            for a in cands:
                if w[a, cur].sum() == 0:
                    break                       # no shared UE: adding cannot reduce the cut
                if zone_2q_cost(inst, cur + [a]) <= budget:
                    cur.append(a)
                    unassigned.discard(a)
                    added = True
                    break
            if not added:
                break
        zones.append(sorted(cur))
        for a in cur:
            zone_of[a] = len(zones) - 1

    # ---- refinement: move a single boundary AP if the cut strictly drops --
    best_cut = _cut(zone_of, ue_aps)
    for _ in range(refine_passes):
        improved = False
        for a in range(n_ap):
            za = zone_of[a]
            if len(zones[za]) <= 1:
                continue
            neigh = {zone_of[b] for b in np.where(w[a] > 0)[0]} - {za}
            for zb in neigh:
                trial = zone_of.copy()
                trial[a] = zb
                if zone_2q_cost(inst, [x for x in range(n_ap) if trial[x] == zb]) > budget:
                    continue
                c = _cut(trial, ue_aps)
                if c < best_cut:
                    zone_of = trial
                    best_cut = c
                    zones[za] = [x for x in zones[za] if x != a]
                    zones[zb] = sorted(zones[zb] + [a])
                    improved = True
                    break
        if not improved:
            break

    # drop empty zones and renumber
    zones = [z for z in zones if z]
    zone_of = np.full(n_ap, -1, dtype=int)
    for zi, z in enumerate(zones):
        for a in z:
            zone_of[a] = zi

    ue_zones = [sorted({int(zone_of[a]) for a in aps}) for aps in ue_aps]
    boundary = np.array([i for i, zs in enumerate(ue_zones) if len(zs) > 1], dtype=int)
    stats = [zone_2q_cost(inst, z, detail=True) for z in zones]

    p = Partition(zone_of=zone_of, zones=zones, ue_zones=ue_zones,
                  boundary=boundary, stats=stats, budget=budget)
    if verbose:
        print(p.summary())
    return p


if __name__ == "__main__":
    from model_instance import make_instance

    for g in (4, 5, 6, 7):
        inst = make_instance(g=g, seed=1, max_deg=2, n_rb_per_ap=4)
        p = partition_aps(inst)
        s = p.summary()
        print(f"g={g}  APs={inst.n_ap:3d} UEs={inst.n_ue:3d} -> zones={s['n_zones']:3d}  "
              f"max2q={s['max_2q']:5d} (budget {p.budget})  maxUE={s['max_ue']:2d}  "
              f"maxQ={s['max_qubits']:2d}  |B|={s['n_boundary']:3d} "
              f"({s['boundary_density']*100:.0f}%)")
