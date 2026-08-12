"""
Centralized strict optimum -- the denominator of the scaling figure.
====================================================================
Solves the joint AP-RB problem of Eq. (joint-problem) exactly:

    max  sum_i sum_{r in V_i} u[i,r] y[i,r]
    s.t. sum_{r in V_i} y[i,r] = 1          for every UE      (one choice)
         sum_i y[i,r]         <= W_r        for every RB      (RB access)
         sum_{i, r in V_i cap R_a} y[i,r] <= W_a   for every AP (AP admission)
         y binary

This is a reference value, not a competing protocol: the distributed result
is reported as a fraction of it.  HiGHS through `scipy.optimize.milp` solves
the sizes used here to proven optimality; the LP relaxation is available as
a cheap valid upper bound if a larger instance is ever needed.
"""

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix


def _matrices(inst):
    """Column layout: one variable per (UE, candidate) pair."""
    cols, obj, ue_of, rb_of = [], [], [], []
    for i, (c, u) in enumerate(zip(inst.cand, inst.util)):
        for r, uu in zip(c, u):
            cols.append((i, int(r)))
            obj.append(float(uu))
            ue_of.append(i)
            rb_of.append(int(r))
    return np.array(ue_of), np.array(rb_of), np.array(obj), cols


def solve_centralized(inst, time_limit=120.0):
    """Return (best_utility, assignment dict UE->RB, status_ok)."""
    ue_of, rb_of, obj, cols = _matrices(inst)
    n = len(obj)

    rows, cidx, vals, lo, hi = [], [], [], [], []
    nr = 0

    # one choice per UE (equality)
    for i in range(inst.n_ue):
        k = np.where(ue_of == i)[0]
        rows += [nr] * len(k); cidx += list(k); vals += [1.0] * len(k)
        lo.append(1.0); hi.append(1.0); nr += 1

    # RB access limit
    for r in np.unique(rb_of):
        k = np.where(rb_of == r)[0]
        rows += [nr] * len(k); cidx += list(k); vals += [1.0] * len(k)
        lo.append(-np.inf); hi.append(float(inst.w_rb)); nr += 1

    # AP admission limit
    ap_of = inst.rb_owner[rb_of]
    for a in np.unique(ap_of):
        k = np.where(ap_of == a)[0]
        rows += [nr] * len(k); cidx += list(k); vals += [1.0] * len(k)
        lo.append(-np.inf); hi.append(float(inst.w_ap)); nr += 1

    A = coo_matrix((vals, (rows, cidx)), shape=(nr, n))
    res = milp(c=-obj, constraints=LinearConstraint(A, lo, hi),
               integrality=np.ones(n), bounds=Bounds(0, 1),
               options=dict(time_limit=time_limit, presolve=True))
    if not res.success:
        return None, None, False
    x = np.round(res.x).astype(int)
    assign = {int(cols[k][0]): int(cols[k][1]) for k in np.where(x == 1)[0]}
    return float(obj @ x), assign, True


def check_feasible(inst):
    """Global feasibility of the instance (needed for every zone to be feasible)."""
    val, _, ok = solve_centralized(inst)
    return ok


if __name__ == "__main__":
    from haiq_instance import make_instance

    for g in (4, 5, 6, 7, 8):
        inst = make_instance(g=g, seed=1, max_deg=2, n_rb_per_ap=4, radius=1.0)
        val, assign, ok = solve_centralized(inst)
        print(f"g={g}  UEs={inst.n_ue:3d}  feasible={ok}  J*={val if val is None else round(val,2)}")
