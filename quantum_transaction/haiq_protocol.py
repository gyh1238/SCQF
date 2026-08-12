"""
Boundary coordination: belief product + confidence-ordered decimation.
======================================================================
Implements Sec. IV-C of the manuscript on top of the zone-local sampler.

Each zone reports one distribution-valued message: its K_z retained joint
draws restricted to its boundary UEs, plus the local utility of each draw.
Coordination is then entirely classical:

    pi_{z,i}(v)   boundary marginal of UE i held by zone z
    b_i(v)        proportional to prod_z pi_{z,i}(v)        (Eq. belief)
    commit        the least ambiguous UE first, then condition every zone's
                  retained list on the committed value and repeat

A commitment is only made to a value that leaves every zone containing that
UE with a non-empty conditional feasible set, so the protocol always
terminates on a globally feasible assignment: by Eq. (factorization), local
strict feasibility plus boundary agreement is global feasibility.  When
conditioning drives a zone's retained list below K_min the zone is
re-sampled with the current commitments pinned -- the exception stage of
Sec. V-D, counted and reported.
"""

import numpy as np

from haiq_zone import (build_zone, enumerate_zone, sample_zone,
                       choose_execution_exponent)


def _marginal(codes, w, k, n_val, alpha=0.5):
    """Boundary marginal of local UE k, weighted and Laplace-smoothed.

    Smoothing matters because a zero is a finite-sample statement, not a
    statement that the value is infeasible; a hard zero would veto a value
    that another zone strongly supports.  The weights are the reconstruction
    weights of Sec. IV-C: they are all equal when the zone executed at the
    target exponent, and unequal when it had to back off.
    """
    c = np.bincount(codes[:, k], weights=w, minlength=n_val).astype(float)
    return (c + alpha) / (c.sum() + alpha * n_val)


def _ess(w):
    """Effective sample size of a weighted list: (sum w)^2 / sum w^2."""
    if len(w) == 0:
        return 0.0
    s1 = float(w.sum())
    s2 = float((w ** 2).sum())
    return s1 * s1 / s2 if s2 > 0 else 0.0


def _weights(util, d_lambda, n_target):
    """Reconstruction weights exp[(lambda - lambda_z) J_z], scaled so that
    they sum to the sample count; that keeps the Laplace smoothing above on
    the same scale whether or not the zone backed off."""
    if len(util) == 0:
        return np.zeros(0)
    lw = d_lambda * np.asarray(util, dtype=float)
    w = np.exp(lw - lw.max())
    s = w.sum()
    return w * (n_target / s) if s > 0 else np.ones(len(w))


def run_protocol(inst, part, beta, ubar, k_accept=200, k_min=25,
                 shot_budget=10_000, rng=None, collect=False, max_attempts=4):
    """
    Run the local stage and the coordination stage.

    Confidence-ordered decimation commits irreversibly, so on rare instances
    it can strand a zone.  The manuscript's own remedy for exhaustion is to
    re-sample, and that is what an attempt is here: every zone is drawn afresh,
    which reorders the beliefs and therefore the commitments.  A run that never
    succeeds is reported with `feasible=False` rather than repaired by some
    other algorithm.
    """
    rng = rng or np.random.default_rng(0)
    last = None
    for attempt in range(max_attempts):
        res = _run_once(inst, part, beta, ubar, k_accept, k_min,
                        shot_budget, rng, collect)
        res["attempts"] = attempt + 1
        last = res
        if res["feasible"]:
            return res
    return last


def _run_once(inst, part, beta, ubar, k_accept, k_min, shot_budget, rng,
              collect):
    """One local stage followed by one coordination pass."""
    lam = beta / ubar

    # ---------------- local stage: one sampling stage per zone -----------
    # Each zone executes at the largest exponent its own shot budget allows
    # and records J_z per draw; the target exponent is restored afterwards by
    # reweighting, so no zone is forced to run at an exponent it cannot pay
    # for and no zone dictates the exponent used for the merge.
    zones, reports = [], []
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue == 0:
            continue
        rows, pref, uref = enumerate_zone(zone, lam)         # exact target law
        beta_z, mu_z, shots, k_eff = choose_execution_exponent(
            zone, beta, ubar, k_accept, shot_budget)
        lam_z = beta_z / ubar
        codes, util, mu_hat, draws = sample_zone(zone, lam_z, k_eff, rng)
        w = _weights(util, lam - lam_z, max(len(util), 1))
        zones.append(zone)
        reports.append(dict(zone=zone, rows=rows, pref=pref, uref=uref,
                            beta_z=beta_z, mu=mu_z, mu_hat=mu_hat,
                            codes=codes, util=util, weights=w,
                            k_eff=k_eff, shots=shots, ess=_ess(w)))

    # ---------------- boundary bookkeeping -------------------------------
    # A UE keeps the same full candidate domain in every zone (full-domain
    # convention), so code indices are directly comparable across zones.
    holders = {}                       # global UE -> [(report index, local k)]
    for ri, rep in enumerate(reports):
        for k in rep["zone"].boundary_local:
            holders.setdefault(rep["zone"].ue[k], []).append((ri, k))
    bnd = [i for i, hs in holders.items() if len(hs) > 1]

    base = [rep["codes"].copy() for rep in reports]   # retained draws per zone
    base_w = [rep["weights"].copy() for rep in reports]   # their target-exponent weights
    base_u = [rep["util"].copy() for rep in reports]      # J_z of each draw
    feas = [rep["rows"] for rep in reports]           # exact F_z per zone

    committed = {}
    exceptions = 0
    backtracks = 0
    trace = []

    def zone_pins(ri, extra=None):
        """Local pin map of a zone: committed values, overridden by `extra`
        (a tentative global UE -> value map used while testing a commitment)."""
        zone = reports[ri]["zone"]
        extra = extra or {}
        pins = {}
        for kk in range(zone.n_ue):
            u = zone.ue[kk]
            if u in extra:
                pins[kk] = int(extra[u])
            elif u in committed:
                pins[kk] = committed[u]
        return pins

    def cond_ok(ri, pins):
        """Is the conditional feasible set of report `ri` non-empty?"""
        rows = feas[ri]
        if len(rows) == 0:
            return False
        m = np.ones(len(rows), dtype=bool)
        for k, v in pins.items():
            m &= rows[:, k] == v
        return bool(m.any())

    def conditioned(ri):
        """Retained draws of zone `ri` consistent with the current commitments,
        with their reconstruction weights and local utilities."""
        zone = reports[ri]["zone"]
        rows = base[ri]
        if len(rows) == 0:
            return rows, base_w[ri], base_u[ri]
        m = np.ones(len(rows), dtype=bool)
        for kk in range(zone.n_ue):
            if zone.ue[kk] in committed:
                m &= rows[:, kk] == committed[zone.ue[kk]]
        return rows[m], base_w[ri][m], base_u[ri][m]

    bnd_set = set(bnd)

    def admissible(i, v, extra=None):
        """Would committing UE i to v leave every zone holding it completable?"""
        e = dict(extra or {})
        e[i] = int(v)
        return all(cond_ok(ri, zone_pins(ri, e)) for ri, k in holders[i])

    def neighbours(i):
        """Uncommitted boundary UEs that share a zone with UE i."""
        s = set()
        for ri, _ in holders[i]:
            zone = reports[ri]["zone"]
            s |= {zone.ue[kk] for kk in range(zone.n_ue)
                  if zone.ue[kk] in bnd_set and zone.ue[kk] not in committed}
        s.discard(i)
        return s

    def safe(i, v):
        """Admissible, and leaving every neighbouring boundary UE a value.

        Checking one step ahead is what keeps confidence-ordered decimation
        from walking into a dead end that only irreversible commitments can
        create; without it the procedure occasionally strands a zone.
        """
        if not admissible(i, v):
            return False
        e = {i: int(v)}
        for j in neighbours(i):
            if not any(admissible(j, w, e) for w in range(len(inst.cand[j]))):
                return False
        return True

    def value_order(i, b):
        """Belief order, with neighbour-preserving values brought to the front."""
        vals = [int(v) for v in np.argsort(-b)]
        good = [v for v in vals if safe(i, v)]
        return good + [v for v in vals if v not in good]

    # Confidence-ordered decimation.  Variables are opened in order of the
    # confidence of their belief, and each is tried against the values its
    # belief prefers -- the greedy path of Sec. IV-C.  Because commitments are
    # irreversible in that description they can strand a later UE, so the
    # search keeps the ability to withdraw the most recent commitment.  With
    # the neighbour check ordering the values, that withdrawal is almost never
    # exercised; `backtracks` reports how often it was.
    seq = []                          # variables in the order they were opened
    order_of, pos_of, belief_of = {}, {}, {}
    max_backtracks = 8 * max(len(bnd), 1)
    aborted = False
    level = 0

    while level < len(bnd):
        if backtracks > max_backtracks:
            aborted = True
            break

        if level == len(seq):
            # open the next variable: the least ambiguous uncommitted UE
            kept = [conditioned(ri) for ri in range(len(reports))]
            best_i, best_conf, best_b = None, -1.0, None
            for i in bnd:
                if i in committed:
                    continue
                n_val = len(inst.cand[i])
                b = np.ones(n_val)
                for ri, k in holders[i]:
                    ck, wk, _ = kept[ri]
                    if len(ck) == 0:
                        continue
                    b *= _marginal(ck, wk, k, n_val)
                b = b / b.sum()
                conf = float(b.max())
                if conf > best_conf:
                    best_i, best_conf, best_b = i, conf, b
            seq.append(best_i)
            order_of[best_i] = value_order(best_i, best_b)
            pos_of[best_i] = -1
            belief_of[best_i] = (best_b.copy(), best_conf)

        i = seq[level]
        p = pos_of[i] + 1
        while p < len(order_of[i]) and not admissible(i, order_of[i][p]):
            p += 1

        if p < len(order_of[i]):
            committed[i] = int(order_of[i][p])
            pos_of[i] = p
            level += 1
            # exception stage: a zone whose conditioned list ran short is re-drawn
            for ri, _ in holders[i]:
                ck, wk, _ = conditioned(ri)
                if _ess(wk) < k_min:
                    rep = reports[ri]
                    zone = rep["zone"]
                    lam_z = rep["beta_z"] / ubar
                    new, new_u, _, _ = sample_zone(zone, lam_z, rep["k_eff"],
                                                   rng, pinned=zone_pins(ri))
                    if len(new):
                        base[ri] = new
                        base_u[ri] = new_u
                        base_w[ri] = _weights(new_u, lam - lam_z, max(len(new_u), 1))
                        exceptions += 1
        else:
            # this variable is exhausted: withdraw it and re-try the previous
            committed.pop(i, None)
            pos_of[i] = -1
            del seq[level:]
            level -= 1
            backtracks += 1
            if level < 0:
                aborted = True
                break
            committed.pop(seq[level], None)

    for u in seq:
        if u in committed:
            b, c = belief_of[u]
            trace.append(dict(ue=u, belief=b, value=committed[u], conf=c))

    # ---------------- final assembly -------------------------------------
    # A zone resolves its interior from the draws it actually holds, not from
    # an enumeration of F_z: with a shot budget it never sees the whole
    # feasible set, and pretending otherwise would hide the cost of the
    # budget.  A zone left with nothing consistent re-runs once with the
    # boundary pinned, which is the same exception stage used above.
    assign = {}
    final_resamples = 0
    for ri, rep in enumerate(reports):
        zone = rep["zone"]
        cand_rows, _, cand_u = conditioned(ri)
        if len(cand_rows) == 0:
            lam_z = rep["beta_z"] / ubar
            cand_rows, cand_u, _, _ = sample_zone(zone, lam_z, rep["k_eff"],
                                                  rng, pinned=zone_pins(ri))
            final_resamples += 1
        if len(cand_rows) == 0:
            continue
        pick = cand_rows[int(np.argmax(cand_u))]
        for k in range(zone.n_ue):
            i = zone.ue[k]
            if i in committed:
                assign[i] = committed[i]
            elif i not in assign:
                assign[i] = int(pick[k])

    # ---------------- scoring --------------------------------------------
    util_total = 0.0
    rb_load, ap_load = {}, {}
    ok_global = True
    for i, ci in assign.items():
        r = int(inst.cand[i][ci])
        util_total += float(inst.util[i][ci])
        rb_load[r] = rb_load.get(r, 0) + 1
        a = int(inst.rb_owner[r])
        ap_load[a] = ap_load.get(a, 0) + 1
    if any(v > inst.w_rb for v in rb_load.values()):
        ok_global = False
    if any(v > inst.w_ap for v in ap_load.values()):
        ok_global = False
    if len(assign) != inst.n_ue:
        ok_global = False

    # ---------------- protocol resources ---------------------------------
    bits = 0
    for rep in reports:
        zone = rep["zone"]
        w = sum(max(1, int(np.ceil(np.log2(len(zone.cand[k])))))
                for k in zone.boundary_local)
        bits += rep["k_eff"] * (w + 32)      # codes + one utility scalar per draw

    out = dict(utility=util_total, feasible=ok_global, assign=assign,
               n_boundary=len(bnd), exceptions=exceptions,
               final_resamples=final_resamples, backtracks=backtracks,
               aborted=aborted, comm_bits=bits,
               shots_max=max((r["shots"] for r in reports), default=0.0),
               shots_total=sum(r["shots"] for r in reports),
               beta_min=min((r["beta_z"] for r in reports), default=beta),
               n_backed_off=sum(1 for r in reports if r["beta_z"] < beta - 1e-6),
               ess_min=min((r["ess"] for r in reports), default=0.0),
               mu_min=min((r["mu"] for r in reports), default=1.0),
               n_active_zones=len(reports))
    if collect:
        out["reports"] = reports
        out["trace"] = trace
        out["holders"] = holders
        out["committed"] = committed
    return out


if __name__ == "__main__":
    from haiq_instance import make_instance, utility_scale
    from haiq_partition import partition_aps
    from haiq_reference import solve_centralized

    for g in (4, 5, 6):
        inst = make_instance(g=g, seed=2)
        opt, _, ok = solve_centralized(inst)
        if not ok:
            print(f"g={g}: infeasible instance, skipped")
            continue
        part = partition_aps(inst)
        res = run_protocol(inst, part, beta=1.5, ubar=utility_scale(inst),
                           rng=np.random.default_rng(1))
        print(f"g={g} UEs={inst.n_ue:3d} zones={res['n_active_zones']:2d} "
              f"|B|={res['n_boundary']:3d}  J={res['utility']:.1f}/{opt:.1f} "
              f"= {100*res['utility']/opt:.1f}%  feasible={res['feasible']}  "
              f"shots<={res['shots_max']:.0f}  backed off={res['n_backed_off']} "
              f"(beta>={res['beta_min']:.2f})  ESS>={res['ess_min']:.0f}")
