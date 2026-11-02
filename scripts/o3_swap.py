"""O3: the drift-swap counterfactual, run to the rule pre-registered in B_SWAP_DESIGN.md.

Both halves are swapped SEPARATELY, so the P2 interleaved-day decoupling survives into the
counterfactual world: the synthetic matcher half derives from u's and v's M-half records only,
the synthetic drift half from their D-half records only. No record crosses halves.

Per user u and donor v:
  synthetic TM  = swap(u <- v) on the M half   -> what the matcher sees
  synthetic TD  = swap(u <- v) on the D half   -> what drift is measured on
  d_app   = JS(real RD, synthetic TD)   the drift actually present in the swapped world
  d_own   = JS(real RD, real TD)        u's own real drift
  d_donor = JS(v's RD, v's TD)          the donor's own drift

Decision rule (B_SWAP_DESIGN §7, fixed before this was run):
  SUPPORTED    the calibration gap tracks d_app and not d_own
  FAILS        it tracks d_own irrespective of d_app -> drift is a marker, not a cause
  INCONCLUSIVE admissible fraction < 50%, or S2 magnitude delivery < 60%, or S1 and S2
               disagree in direction

Usage: o3_swap.py <windows.tsv> [n_users] [seed]
"""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M

PROF = sys.argv[1]
NU   = int(sys.argv[2]) if len(sys.argv) > 2 else 0
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20261102
CAP, GRID = 50.0, 200

U = M.parse_windows(PROF, keep={"W3M", "W3D", "W4M", "W4D"})
REN = {"W3M": "RM", "W3D": "RD", "W4M": "TM", "W4D": "TD"}
U = {u: {REN[g]: v for g, v in gs.items()} for u, gs in U.items()}
users = sorted([u for u, g in U.items() if len(g) == 4], key=int)
rng = np.random.default_rng(SEED)
if NU and NU < len(users):
    users = [users[i] for i in sorted(rng.choice(len(users), NU, replace=False))]
N = len(users)
UIDX = {u: i for i, u in enumerate(users)}
print(f"Arm A, users with all four (window x half) groups: {N:,}\n")

# ---------------------------------------------------------------- offset frames, per half
def xy(k): a, b = k.split("_"); return int(a), int(b)
def anchor(c): return xy(max(c.items(), key=lambda kv: (kv[1], kv[0]))[0])
def offs(c, a):
    tot = sum(c.values())
    out = {}
    for k, v in c.items():
        x, y = xy(k); out[(x - a[0], y - a[1])] = v/tot
    return out

P = {}
same_anchor = 0
for u in users:
    e = {}
    for half in "MD":
        aR = anchor(U[u]["R"+half]["c"]); aT = anchor(U[u]["T"+half]["c"])
        e[half] = dict(aR=aR, D=(aT[0]-aR[0], aT[1]-aR[1]),
                       qR=offs(U[u]["R"+half]["c"], aR),
                       qT=offs(U[u]["T"+half]["c"], aT))
    same_anchor += int(e["M"]["aR"] == e["D"]["aR"])
    P[u] = e
print(f"the M-half and D-half reference anchors agree for {100*same_anchor/N:.1f}% of users "
      "(they are estimated separately to keep the halves disjoint)\n")

# population drift field, per half, for the POP condition
POP = {}
for half in "MD":
    QR, QT = {}, {}
    for u in users:
        for o, m in P[u][half]["qR"].items(): QR[o] = QR.get(o, 0.0) + m
        for o, m in P[u][half]["qT"].items(): QT[o] = QT.get(o, 0.0) + m
    sR = sum(QR.values()); sT = sum(QT.values())
    QR = {o: m/sR for o, m in QR.items()}; QT = {o: m/sT for o, m in QT.items()}
    nov = sum(m for o, m in QT.items() if o not in QR)
    POP[half] = dict(qR=QR, qT=QT, m_nov=nov)
print(f"population drift field: novelty mass {POP['M']['m_nov']:.3f} (M half), "
      f"{POP['D']['m_nov']:.3f} (D half)\n")

def jsd(p, q):
    keys = sorted(set(p) | set(q))
    a = np.array([p.get(k, 0.0) for k in keys]); b = np.array([q.get(k, 0.0) for k in keys])
    return M.js(a, b)

def swap(u, half, donor=None, novelty="radius", pop=False):
    """Return (synthetic offset distribution, synthetic anchor) or (None, None)."""
    pu = P[u][half]; qR = pu["qR"]
    if pop:
        qvR, qvT, m_nov = POP[half]["qR"], POP[half]["qT"], POP[half]["m_nov"]
        Dv = (0, 0)
    else:
        pv = P[donor][half]; qvR, qvT = pv["qR"], pv["qT"]
        m_nov = sum(m for o, m in qvT.items() if o not in qvR); Dv = pv["D"]
    nov = {o: m for o, m in qvT.items() if o not in qvR}
    shared = [o for o in qR if o in qvR]
    if not shared: return None, None
    ratio = {o: min(CAP, qvT.get(o, 0.0)/qvR[o]) for o in shared}
    kappa = float(np.mean(list(ratio.values())))
    w = {o: qR[o]*ratio.get(o, kappa) for o in qR}
    sw = sum(w.values())
    if sw <= 0: return None, None
    out = {o: (1.0-m_nov)*x/sw for o, x in w.items()}
    if nov and m_nov > 0:
        if novelty == "exact":
            t = sum(nov.values())
            for o, m in nov.items(): out[o] = out.get(o, 0.0) + m_nov*m/t
        else:
            radii = np.array([math.hypot(*o) for o in nov])
            mass = np.array(list(nov.values())); mass = mass/mass.sum()
            pick = rng.choice(len(radii), len(nov), replace=True, p=mass)
            taken = set(out); drawn = {}
            for j in pick:
                rr = radii[j]
                for _ in range(16):
                    th = rng.uniform(0, 2*math.pi)
                    o = (int(round(rr*math.cos(th))), int(round(rr*math.sin(th))))
                    if o not in taken and o not in drawn: drawn[o] = 1.0; break
                else: drawn[o] = drawn.get(o, 0.0) + 1.0
            t = sum(drawn.values()) or 1.0
            for o, m in drawn.items(): out[o] = out.get(o, 0.0) + m_nov*m/t
    out = {o: m for o, m in out.items() if m > 1e-12}
    if not out: return None, None
    return out, (pu["aR"][0]+Dv[0], pu["aR"][1]+Dv[1])

POPMAG = {}
def swap_S1(u, half, donor):
    """Magnitude-only: interpolate u's own ref->tgt to a target drift magnitude -- the
    donor's own drift, or the population median when donor is None (the POP rung)."""
    pu = P[u][half]
    target = POPMAG[half] if donor is None else jsd(P[donor][half]["qR"], P[donor][half]["qT"])
    keys = sorted(set(pu["qR"]) | set(pu["qT"]))
    lo, hi = 0.0, 1.0
    for _ in range(30):
        mid = 0.5*(lo+hi)
        mix = {k: (1-mid)*pu["qR"].get(k, 0.0) + mid*pu["qT"].get(k, 0.0) for k in keys}
        if jsd(pu["qR"], mix) < target: lo = mid
        else: hi = mid
    mix = {k: v for k, v in
           ((k, (1-hi)*pu["qR"].get(k, 0.0) + hi*pu["qT"].get(k, 0.0)) for k in keys)
           if v > 1e-12}
    aT = (pu["aR"][0]+P[u][half]["D"][0], pu["aR"][1]+P[u][half]["D"][1])
    return mix, aT

def to_cells(q, a, total):
    """Offsets -> absolute cell counts, dropping mass that leaves the 200x200 grid."""
    c = {}; drop = 0.0
    for (dx, dy), m in q.items():
        x, y = a[0]+dx, a[1]+dy
        if 0 <= x < GRID and 0 <= y < GRID: c[f"{x}_{y}"] = c.get(f"{x}_{y}", 0.0) + m
        else: drop += m
    s = sum(c.values())
    if s <= 0: return None, 1.0
    return {k: v*total/s for k, v in c.items()}, drop

# ---------------------------------------------------------------- donors
base = np.array([[len(P[u]["D"]["qR"]),
                  M.ent(np.array(list(P[u]["D"]["qR"].values()))),
                  math.log(U[u]["RD"]["n"] + U[u]["TD"]["n"])] for u in users])
bz = (base - base.mean(0))/base.std(0)
order = np.lexsort((bz[:, 2], bz[:, 1], bz[:, 0]))
nn = np.empty(N, int)
for r, i in enumerate(order):
    nn[i] = order[r+1] if r+1 < N else order[r-1]
don_nn = [users[j] for j in nn]
don_rnd = [users[j] for j in rng.permutation(N)]
print(f"donors: NN matched on (cell count, offset entropy, log coverage) via lexsorted "
      f"neighbour; mean |z| distance {np.abs(bz - bz[nn]).mean():.3f}")
print(f"        RND is a random permutation; self-donor collisions in RND: "
      f"{sum(1 for i in range(N) if don_rnd[i]==users[i])}\n")

# ---------------------------------------------------------------- plausibility gate
def marg(q):
    w = np.array(list(q.values())); w = w/w.sum()
    r = np.array([math.hypot(*o) for o in q])
    return (len(q), M.ent(w), float(w.max()), float((r*w).sum()))
covD = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)
dec = np.digitize(covD, np.quantile(covD, np.arange(.1, 1.0, .1)))
real = np.array([marg(P[u]["D"]["qT"]) for u in users], float)
BAND = {}
for k in range(10):
    m = dec == k
    BAND[k] = (np.quantile(real[m], .025, axis=0), np.quantile(real[m], .975, axis=0))
print("plausibility gate: the four marginals must lie in the real target 95% band of the")
print("user's OWN coverage decile (B_SWAP_DESIGN §5 required the coverage-matched version)\n")

# ---------------------------------------------------------------- run the ladder
cal = rng.permutation(N); ev = cal[N//2:]; cal = cal[:N//2]
def auc(p, y):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1*(n1+1)/2)/(n1*n0)) if n1 and n0 else float("nan")
def z(v): return (v - v.mean())/max(v.std(), 1e-12)

REALCAL = {}
d_own = np.array([jsd(P[u]["D"]["qR"], P[u]["D"]["qT"]) for u in users])
for half in "MD":
    POPMAG[half] = float(np.median([jsd(P[u][half]["qR"], P[u][half]["qT"]) for u in users]))
print(f"population median drift magnitude: M {POPMAG['M']:.3f}, D {POPMAG['D']:.3f}\n")
RESULT = {}
CONDS = [("REAL  (no swap)",        None,    None,  False, None),
         ("S3    self, exact",      "self",  "exact", False, "S2"),
         ("S2    NN donor",         don_nn,  "radius", False, "S2"),
         ("S2    RND donor",        don_rnd, "radius", False, "S2"),
         ("S2    POP field",        None,    "radius", True,  "S2"),
         ("S1    NN donor",         don_nn,  None,   False, "S1"),
         ("S1    POP magnitude",     "popmag", None,   False, "S1")]

for label, donors, novmode, pop, kind in CONDS:
    Uc = {}; d_app = np.full(N, np.nan); d_don = np.full(N, np.nan)
    adm = np.zeros(N, bool); dropped = []
    for i, u in enumerate(users):
        g = {"RM": U[u]["RM"], "RD": U[u]["RD"]}
        ok = True
        syn = {}
        for half in "MD":
            if kind is None:
                q, a = P[u][half]["qT"], (P[u][half]["aR"][0]+P[u][half]["D"][0],
                                          P[u][half]["aR"][1]+P[u][half]["D"][1])
            elif kind == "S1":
                q, a = swap_S1(u, half, None if donors == "popmag" else donors[i])
            else:
                dv = u if donors == "self" else (None if pop else donors[i])
                q, a = swap(u, half, donor=dv, novelty=novmode, pop=pop)
            if q is None: ok = False; break
            syn[half] = (q, a)
        if not ok: continue
        for half in "MD":
            q, a = syn[half]
            c, dr = to_cells(q, a, U[u]["T"+half]["n"])
            if c is None: ok = False; break
            dropped.append(dr)
            g["T"+half] = dict(n=U[u]["T"+half]["n"], d=U[u]["T"+half]["d"], n27=0,
                               c=c, s=U[u]["T"+half]["s"])
        if not ok: continue
        Uc[u] = g
        d_app[i] = jsd(P[u]["D"]["qR"], syn["D"][0])
        if kind is not None and donors not in (None, "self", "popmag"):
            d_don[i] = jsd(P[donors[i]]["D"]["qR"], P[donors[i]]["D"]["qT"])
        mg = marg(syn["D"][0]); lo, hi = BAND[dec[i]]
        adm[i] = bool(np.all(np.array(mg) >= lo) and np.all(np.array(mg) <= hi))
    sub = [u for u in users if u in Uc]
    idx = np.array([UIDX[u] for u in sub])
    res = M.match(Uc, sub, "RM", "TM")
    cv = covD[idx]
    ii = {u: j for j, u in enumerate(sub)}
    calm = np.array([ii[users[i]] for i in cal if users[i] in ii])
    evm  = np.array([ii[users[i]] for i in ev  if users[i] in ii])
    k = M.calibrate(res, cv, calm, evm)
    y, p = res["y"], k["p_platt"]
    # AMENDMENT 2026-11-02 (A_PREMISES s8): also carry the REAL world's calibrator into this
    # world unchanged. The pre-registered plan refitted the calibrator inside the swapped
    # world, which absorbs the applied drift and leaves nothing for the decision test to
    # measure. The transported variant mirrors O1 s4b, which is where the effect lives.
    if label.startswith("REAL"):
        REALCAL["w"], REALCAL["mu"], REALCAL["sd"] = k["w"], k["mu"], k["sd"]
        p_tr = p
    else:
        X = M.platt_features(res["topk"], cv)
        Xz = (X - REALCAL["mu"]) / REALCAL["sd"]
        p_tr = M.logit_apply(Xz, REALCAL["w"])
    da, do = d_app[idx], d_own[idx]
    am = adm[idx]
    RESULT[label] = dict(sub=sub, idx=idx, y=y, p=p, p_tr=p_tr, da=da, do=do, dd=d_don[idx],
                         adm=am, ev=evm, drop=float(np.mean(dropped)) if dropped else 0.0)
    print(f"{label:<20} n={len(sub):>6,}  top-1 {y.mean():.4f}  ECE(ev) {M.ece(p[evm],y[evm]):.4f}  "
          f"ECE_transported {M.ece(p_tr[evm],y[evm]):.4f}  "
          f"mean d_app {np.nanmean(da):.3f}  admissible {100*am.mean():.1f}%  "
          f"grid-loss {100*np.mean(dropped):.2f}%")
np.savez("results/o3_vars.npz", **{f"{k.split()[0]}_{k.split()[1]}_{f}": v
         for k, d in RESULT.items() for f, v in d.items()
         if isinstance(v, np.ndarray)})
import pickle
pickle.dump({k: {kk: vv for kk, vv in d.items() if kk != "sub"} for k, d in RESULT.items()},
            open("results/o3_result.pkl", "wb"))
print("\nsaved results/o3_vars.npz and results/o3_result.pkl")
