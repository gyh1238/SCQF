"""
Coverage-instance generator for the joint AP-RB assignment model.
=================================================================
Produces the data the manuscript's system model (Sec. III) assumes:

    A            AP set, each AP `a` owning a disjoint RB pool R_a
    U            UE set, each UE `i` with a coverage-dependent candidate
                 set V_i subseteq R  (the RB variable is the decision;
                 the serving AP is the owner of the selected RB)
    u[i, r]      utility of assigning UE i to RB r
    W_r, W_a     RB access limit and AP admission limit

The region is a perturbed AP grid of side `g` (spacing 1.0).  Holding the
AP grid spacing, the UE density, the coverage radius and the candidate
degree fixed while increasing `g` grows the *global* problem without
changing local density -- this is the growth axis used by the scaling
figure, and it is exactly the "bounded zone density and bounded candidate
degree" regime the manuscript invokes in Sec. V-D.
"""

from dataclasses import dataclass
from math import ceil, log2

import numpy as np


@dataclass
class Instance:
    """One coverage instance."""
    ap_xy: np.ndarray          # (n_ap, 2) AP positions
    ue_xy: np.ndarray          # (n_ue, 2) UE positions
    rb_owner: np.ndarray       # (n_rb,) owner AP index of each RB
    cand: list                 # cand[i] = array of candidate RB indices of UE i
    util: list                 # util[i] = array of utilities, aligned with cand[i]
    w_rb: int                  # RB access limit  (per RB)
    w_ap: int                  # AP admission limit (per AP)
    g: int                     # grid side (region size)
    seed: int

    # ---- derived sizes ----
    @property
    def n_ap(self):
        return len(self.ap_xy)

    @property
    def n_ue(self):
        return len(self.ue_xy)

    @property
    def n_rb(self):
        return len(self.rb_owner)

    def code_width(self, i):
        """l_i = ceil(log2 |V_i|)  -- state qubits held for UE i, Eq. (state-width)."""
        return max(1, ceil(log2(len(self.cand[i]))))

    def ues_of_ap(self, a):
        """UEs holding at least one candidate RB owned by AP a."""
        out = []
        for i, c in enumerate(self.cand):
            if np.any(self.rb_owner[c] == a):
                out.append(i)
        return out


def make_instance(g=5, ue_per_ap=2.0, n_rb_per_ap=4, radius=1.2,
                  max_deg=2, w_rb=1, w_ap=4, seed=0, jitter=0.18):
    """
    Build one instance on a g x g perturbed AP grid.

    Densities (UEs per AP, RBs per AP, coverage radius, candidate degree)
    are independent of `g`, so `g` is a pure problem-size knob.
    """
    rng = np.random.default_rng(seed)

    # --- APs on a jittered unit grid -------------------------------------
    gx, gy = np.meshgrid(np.arange(g), np.arange(g))
    ap_xy = np.stack([gx.ravel(), gy.ravel()], axis=1).astype(float)
    ap_xy += rng.normal(0.0, jitter, ap_xy.shape)
    n_ap = len(ap_xy)

    # --- RB pools: AP a owns RBs [a*n_rb_per_ap, (a+1)*n_rb_per_ap) ------
    rb_owner = np.repeat(np.arange(n_ap), n_rb_per_ap)
    n_rb = len(rb_owner)
    # per-RB quality jitter so that RBs of the same AP are not interchangeable
    rb_gain = rng.normal(0.0, 0.6, n_rb)

    # --- UEs uniform over the region -------------------------------------
    n_ue = int(round(ue_per_ap * n_ap))
    ue_xy = rng.uniform(-0.5, g - 0.5, size=(n_ue, 2))

    # --- candidate sets ---------------------------------------------------
    # An AP offers each covered UE one RB from its own pool, allocated round
    # robin by proximity, so that neighbouring UEs are not all pushed onto the
    # same RB.  A UE therefore holds one candidate per covering AP, and it is
    # coverage overlap -- not RB multiplicity inside one AP -- that creates
    # boundary UEs.
    dist = np.linalg.norm(ue_xy[:, None, :] - ap_xy[None, :, :], axis=2)
    offer = {}                                          # offer[(i, a)] = RB index
    for a in range(n_ap):
        covered = np.where(dist[:, a] <= radius)[0]
        covered = covered[np.argsort(dist[covered, a])]
        pool = np.where(rb_owner == a)[0]
        for rank, i in enumerate(covered):
            offer[(int(i), a)] = int(pool[rank % len(pool)])

    cand, util, keep = [], [], []
    for i in range(n_ue):
        near = np.where(dist[i] <= radius)[0]
        if len(near) == 0:
            continue                                   # uncovered UE: drop
        rbs = np.array([offer[(i, int(a))] for a in near])
        u = (10.0 * np.exp(-dist[i, near] / (0.6 * radius)) + rb_gain[rbs]
             + rng.normal(0, 0.15, len(near)))
        u = np.maximum(u, 0.05)
        order = np.argsort(-u)[:max_deg]               # keep the strongest APs
        cand.append(rbs[order])
        util.append(u[order])
        keep.append(i)

    return Instance(ap_xy=ap_xy, ue_xy=ue_xy[keep], rb_owner=rb_owner,
                    cand=cand, util=util, w_rb=w_rb, w_ap=w_ap, g=g, seed=seed)


def utility_scale(inst):
    """
    Global utility scale ubar of Sec. IV-C, used as lambda = beta / ubar.

    It is a property of the utility *model*, not of an instance, so it is
    defined from the model's dynamic range and shared by every zone.
    """
    return float(np.mean([u.max() for u in inst.util]))


if __name__ == "__main__":
    for g in (3, 5, 7):
        inst = make_instance(g=g, seed=1)
        degs = [len(c) for c in inst.cand]
        print(f"g={g:2d}  APs={inst.n_ap:3d}  RBs={inst.n_rb:4d}  UEs={inst.n_ue:3d}  "
              f"deg mean={np.mean(degs):.2f} max={max(degs)}  ubar={utility_scale(inst):.2f}")
