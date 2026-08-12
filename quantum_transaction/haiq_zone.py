"""
Zone-local sampler: the law the accepted branch of the circuit produces.
========================================================================
From Eq. (theta)-(gibbs-local) of the manuscript, the circuit prepares a
uniform superposition over valid codewords, applies one controlled-R_Y per
candidate with

    g_{z,i}(r) = exp[-lambda (u^max_{z,i} - u^{(z)}_{i,r})],
    theta      = 2 arccos sqrt(g),

and accepts on (all cost qubits |0>) AND (all constraint flags 0).  The
accepted-branch law is therefore

    p_z(r_z)  proportional to  q_0(r_z) * prod_i g_{z,i}(r_i) * 1[r_z in F_z]
              proportional to  exp[lambda J_z(r_z)]           on F_z,

which is *exactly* independent per-UE sampling with weight g_{z,i} followed
by rejection on strict feasibility.  The classical routine below is that
same sampler, so it is not an approximation of the quantum output but the
identical law; `haiq_certify.py` verifies the equality against a Qiskit
statevector.  The quantum and classical versions differ only in cost: the
classical sampler needs ~K/mu_z draws for K accepted samples, while
amplitude amplification reaches the accepted branch in O(mu_z^{-1/2})
rounds -- the quadratic factor is the whole of the local quantum gain.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass
class Zone:
    idx: int
    aps: list                  # APs owned
    ue: list                   # U_z, global UE indices
    cand: list                 # cand[k] = candidate RBs of local UE k (full domain)
    uloc: list                 # uloc[k] = owner-local utility, 0 for out-of-zone RBs
    owned: list                # owned[k] = bool mask, candidate owned by this zone
    rbs: list                  # RBs owned by the zone that some local UE can take
    w_rb: int
    w_ap: int
    rb_owner: np.ndarray
    boundary_local: list       # local positions of boundary UEs
    boundary_global: list      # global indices of the same UEs

    @property
    def n_ue(self):
        return len(self.ue)

    def code_widths(self):
        return [max(1, int(np.ceil(np.log2(len(c))))) for c in self.cand]


def build_zone(inst, part, z):
    """Assemble the local view of zone `z` (Sec. III-C, full-domain convention)."""
    aps = part.zones[z]
    ap_set = set(aps)
    in_zone_rb = np.array([o in ap_set for o in inst.rb_owner])

    ue, cand, uloc, owned = [], [], [], []
    for i, c in enumerate(inst.cand):
        m = in_zone_rb[c]
        if not m.any():
            continue
        ue.append(i)
        cand.append(np.asarray(c))
        owned.append(m)
        uloc.append(np.where(m, inst.util[i], 0.0))     # u^{(z)}: 0 outside the zone

    rbs = sorted({int(r) for c, m in zip(cand, owned) for r in c[m]})
    bl = [k for k, i in enumerate(ue) if len(part.ue_zones[i]) > 1]
    return Zone(idx=z, aps=list(aps), ue=ue, cand=cand, uloc=uloc, owned=owned,
                rbs=rbs, w_rb=inst.w_rb, w_ap=inst.w_ap, rb_owner=inst.rb_owner,
                boundary_local=bl, boundary_global=[ue[k] for k in bl])


# ---------------------------------------------------------------- feasibility

def feasible(zone, codes):
    """phi_z = 1: every owned RB within its access limit, every owned AP within
    its admission limit.  Candidates owned elsewhere load neither."""
    rb_load, ap_load = {}, {}
    for k, ci in enumerate(codes):
        if not zone.owned[k][ci]:
            continue                          # value owned outside the zone
        r = int(zone.cand[k][ci])
        rb_load[r] = rb_load.get(r, 0) + 1
        if rb_load[r] > zone.w_rb:
            return False
        a = int(zone.rb_owner[r])
        ap_load[a] = ap_load.get(a, 0) + 1
        if ap_load[a] > zone.w_ap:
            return False
    return True


def local_utility(zone, codes):
    """J_z of Eq. (local-utility)."""
    return float(sum(zone.uloc[k][ci] for k, ci in enumerate(codes)))


# ------------------------------------------------------------ exact zone law

def enumerate_zone(zone, lam, pinned=None):
    """
    Exact enumeration of F_z with the Gibbs weights of Eq. (gibbs-local).

    `pinned` maps local UE position -> code index (committed boundary values).
    Returns (codes array, probabilities, utilities).  Zone widths are bounded
    by the partitioner, so this is a cheap exact reference.
    """
    pinned = pinned or {}
    if zone.n_ue == 0:
        return np.zeros((0, 0), dtype=int), np.zeros(0), np.zeros(0)
    doms = [[pinned[k]] if k in pinned else range(len(zone.cand[k]))
            for k in range(zone.n_ue)]
    rows, util = [], []
    for codes in product(*doms):
        if feasible(zone, codes):
            rows.append(codes)
            util.append(local_utility(zone, codes))
    if not rows:
        return np.zeros((0, zone.n_ue), dtype=int), np.zeros(0), np.zeros(0)
    util = np.asarray(util)
    logp = lam * util
    p = np.exp(logp - logp.max())
    return np.asarray(rows, dtype=int), p / p.sum(), util


def acceptance_mass(zone, lam):
    """
    mu_z of Eq. (mu): the probability that one circuit execution lands in the
    accepted branch, before amplification.  With a uniform preparation over
    valid codewords, mu_z = sum_{F_z} prod_i g_{z,i}(r_i)/|V_i|.
    """
    wts = _per_ue_weights(zone, lam)
    acc = 0.0
    for codes in product(*[range(len(c)) for c in zone.cand]):
        if not feasible(zone, codes):
            continue
        w = 1.0
        for k, ci in enumerate(codes):
            w *= wts[k][ci] / len(zone.cand[k])
        acc += w
    return acc


def _per_ue_weights(zone, lam):
    """g_{z,i}(r) for every candidate of every local UE, Eq. (theta)."""
    out = []
    for k in range(zone.n_ue):
        u = zone.uloc[k]
        out.append(np.exp(-lam * (u.max() - u)))
    return out


# --------------------------------------------------------------- the sampler

def sample_zone(zone, lam, k_accept, rng, pinned=None, max_draws=400_000):
    """
    Draw `k_accept` accepted samples by the rejection process the accepted
    branch implements.  Returns (codes, utilities, mu_hat, n_draws).

    Early abort on the first violated constraint leaves the accepted law
    unchanged and only saves work.
    """
    pinned = pinned or {}
    n = zone.n_ue
    if n == 0:                                         # AP with no candidate UE
        return np.zeros((0, 0), dtype=int), np.zeros(0), 1.0, 0
    wts = _per_ue_weights(zone, lam)
    probs = [w / w.sum() for w in wts]                 # uniform prior x g, normalized

    # per-candidate lookup tables for vectorized constraint checking
    rb_tab = [np.where(zone.owned[k], zone.cand[k], -1) for k in range(n)]
    ap_tab = [np.where(zone.owned[k], zone.rb_owner[zone.cand[k]], -1) for k in range(n)]
    u_tab = [np.asarray(zone.uloc[k]) for k in range(n)]

    kept_codes, kept_util = [], []
    n_kept = 0
    draws = 0
    batch = 4096
    while n_kept < k_accept and draws < max_draws:
        b = min(batch, max_draws - draws)
        if b <= 0:
            break
        codes = np.empty((b, n), dtype=int)
        for k in range(n):
            if k in pinned:
                codes[:, k] = pinned[k]
            else:
                codes[:, k] = rng.choice(len(probs[k]), size=b, p=probs[k])
        draws += b

        rb_sel = np.stack([rb_tab[k][codes[:, k]] for k in range(n)], axis=1)
        ap_sel = np.stack([ap_tab[k][codes[:, k]] for k in range(n)], axis=1)
        ok = np.ones(b, dtype=bool)
        for r in zone.rbs:
            ok &= (rb_sel == r).sum(axis=1) <= zone.w_rb
        for a in zone.aps:
            ok &= (ap_sel == a).sum(axis=1) <= zone.w_ap

        sel = codes[ok]
        n_kept += len(sel)
        if len(sel):
            util = np.stack([u_tab[k][sel[:, k]] for k in range(n)], axis=1).sum(axis=1)
            kept_codes.append(sel)
            kept_util.append(util)

    if kept_codes:
        codes_out = np.concatenate(kept_codes)[:k_accept]
        util_out = np.concatenate(kept_util)[:k_accept]
    else:
        codes_out = np.zeros((0, n), dtype=int)
        util_out = np.zeros(0)
    mu_hat = n_kept / max(draws, 1)
    return codes_out, util_out, mu_hat, draws


def choose_execution_exponent(zone, beta_target, ubar, k_accept, shot_budget,
                              tol=1e-3):
    """
    Largest exponent this zone can afford to execute, Sec. IV-C of the paper.

    A tight zone has a small acceptance mass, and collecting K_z accepted
    draws at the target exponent can cost far more shots than any run has.
    The remedy is not to shrink the target but to execute at a lower
    beta_z <= beta, where the accepted law is flatter and therefore cheaper
    to hit, and to restore the target exponent after measurement by
    reweighting each recorded draw (see `haiq_protocol`).

    Returns (beta_z, mu_z, shots, k_eff).  `k_eff` is below `k_accept` only
    when even the uniform law (beta_z = 0) cannot meet the budget, in which
    case the zone reports fewer draws rather than exceeding it.
    """
    def shots_at(b):
        mu = acceptance_mass(zone, b / ubar)
        return mu, (k_accept / mu if mu > 0 else np.inf)

    mu, s = shots_at(beta_target)
    if s <= shot_budget:
        return beta_target, mu, s, k_accept

    mu0, s0 = shots_at(0.0)
    if s0 > shot_budget:                       # even the uniform law is too dear
        k_eff = max(1, int(shot_budget * mu0))
        return 0.0, mu0, k_eff / mu0, k_eff

    lo, hi = 0.0, beta_target                  # feasible at lo, not at hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if shots_at(mid)[1] <= shot_budget:
            lo = mid
        else:
            hi = mid
    mu, s = shots_at(lo)
    return lo, mu, s, k_accept


def amplification_rounds(mu):
    """k = round(pi/4 / arcsin sqrt(mu)) - 1/2, the standard Grover schedule."""
    mu = float(np.clip(mu, 1e-12, 1.0))
    return max(0, int(round(np.pi / (4 * np.arcsin(np.sqrt(mu))) - 0.5)))


if __name__ == "__main__":
    from haiq_instance import make_instance, utility_scale
    from haiq_partition import partition_aps

    inst = make_instance(g=5, seed=1, max_deg=2, n_rb_per_ap=4)
    part = partition_aps(inst)
    ubar = utility_scale(inst)
    lam = 1.5 / ubar
    rng = np.random.default_rng(0)
    print(f"ubar={ubar:.2f}  lambda={lam:.4f}  zones={part.n_zones}")
    shown = 0
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue == 0 or shown >= 5:
            continue
        shown += 1
        rows, p, u = enumerate_zone(zone, lam)
        codes, ut, mu, dr = sample_zone(zone, lam, 400, rng)
        print(f"  zone {z}: N_z={zone.n_ue} |F_z|={len(rows)} |B_z|={len(zone.boundary_local)} "
              f"mu_hat={mu:.3f} rounds={amplification_rounds(mu)} draws={dr}")
