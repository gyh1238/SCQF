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
                       acceptance_mass, amplification_rounds)


def _marginal(codes, k, n_val, alpha=0.5):
    """Boundary marginal of local UE k with Laplace smoothing.

    Smoothing matters because a zero is a finite-sample statement, not a
    statement that the value is infeasible; a hard zero would veto a value
    that another zone strongly supports.
    """
    c = np.bincount(codes[:, k], minlength=n_val).astype(float)
    return (c + alpha) / (c.sum() + alpha * n_val)


def run_protocol(inst, part, beta, ubar, k_accept=400, k_min=40, rng=None,
                 collect=False, max_attempts=4):
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
        res = _run_once(inst, part, beta, ubar, k_accept, k_min, rng, collect)
        res["attempts"] = attempt + 1
        last = res
        if res["feasible"]:
            return res
    return last


def _run_once(inst, part, beta, ubar, k_accept, k_min, rng, collect):
    """One local stage followed by one coordination pass."""
    lam = beta / ubar

    # ---------------- local stage: one sampling stage per zone -----------
    zones, reports = [], []
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue == 0:
            continue
        rows, pref, uref = enumerate_zone(zone, lam)         # exact zone law
        mu = acceptance_mass(zone, lam)
        k_rounds = amplification_rounds(mu)
        p_acc = np.sin((2 * k_rounds + 1) * np.arcsin(np.sqrt(max(mu, 1e-12)))) ** 2
        codes, util, mu_hat, draws = sample_zone(zone, lam, k_accept, rng)
        zones.append(zone)
        reports.append(dict(zone=zone, rows=rows, pref=pref, uref=uref,
                            mu=mu, mu_hat=mu_hat, rounds=k_rounds,
                            p_acc=p_acc, codes=codes, util=util,
                            shots_quantum=k_accept / max(p_acc, 1e-12),
                            draws_classical=k_accept / max(mu, 1e-12)))

    # ---------------- boundary bookkeeping -------------------------------
    # A UE keeps the same full candidate domain in every zone (full-domain
    # convention), so code indices are directly comparable across zones.
    holders = {}                       # global UE -> [(report index, local k)]
    for ri, rep in enumerate(reports):
        for k in rep["zone"].boundary_local:
            holders.setdefault(rep["zone"].ue[k], []).append((ri, k))
    bnd = [i for i, hs in holders.items() if len(hs) > 1]

    base = [rep["codes"].copy() for rep in reports]   # retained draws per zone
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
        """Retained draws of zone `ri` consistent with the current commitments."""
        zone = reports[ri]["zone"]
        rows = base[ri]
        if len(rows) == 0:
            return rows
        m = np.ones(len(rows), dtype=bool)
        for kk in range(zone.n_ue):
            if zone.ue[kk] in committed:
                m &= rows[:, kk] == committed[zone.ue[kk]]
        return rows[m]

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
                    if len(kept[ri]) == 0:
                        continue
                    b *= _marginal(kept[ri], k, n_val)
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
                if len(conditioned(ri)) < k_min:
                    zone = reports[ri]["zone"]
                    new, _, _, _ = sample_zone(zone, lam, k_accept, rng,
                                               pinned=zone_pins(ri))
                    if len(new):
                        base[ri] = new
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
    assign = {}
    for ri, rep in enumerate(reports):
        zone = rep["zone"]
        rows = feas[ri]
        m = np.ones(len(rows), dtype=bool)
        for k in range(zone.n_ue):
            if zone.ue[k] in committed:
                m &= rows[:, k] == committed[zone.ue[k]]
        cand_rows = rows[m]
        if len(cand_rows) == 0:
            continue
        util = np.array([sum(zone.uloc[k][c[k]] for k in range(zone.n_ue))
                         for c in cand_rows])
        pick = cand_rows[int(np.argmax(util))]
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
        bits += k_accept * (w + 32)          # codes + one utility scalar per draw
    shots_q = sum(r["shots_quantum"] for r in reports)
    draws_c = sum(r["draws_classical"] for r in reports)

    out = dict(utility=util_total, feasible=ok_global, assign=assign,
               n_boundary=len(bnd), exceptions=exceptions, backtracks=backtracks,
               aborted=aborted,
               comm_bits=bits, shots_quantum=shots_q, draws_classical=draws_c,
               max_rounds=max((r["rounds"] for r in reports), default=0),
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
              f"exc={res['exceptions']}  maxk={res['max_rounds']}")
