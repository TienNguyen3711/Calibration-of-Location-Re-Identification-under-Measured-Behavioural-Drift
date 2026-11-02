import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M

PROF = sys.argv[1]
NU   = int(sys.argv[2]) if len(sys.argv) > 2 else 0
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20261102
CAP, GRID = 50.0, 200
TARGET_TOP1 = 0.4943      # S2-NN, from results/O3_swap.md
S2_ECE_TR   = 0.1609      # S2-NN transported-calibrator ECE, same source

U = M.parse_windows(PROF, keep={"W3M", "W3D", "W4M", "W4D"})
REN = {"W3M": "RM", "W3D": "RD", "W4M": "TM", "W4D": "TD"}
U = {u: {REN[g]: v for g, v in gs.items()} for u, gs in U.items()}
users = sorted([u for u, g in U.items() if len(g) == 4], key=int)
rng = np.random.default_rng(SEED)
if NU and NU < len(users):
    users = [users[i] for i in sorted(rng.choice(len(users), NU, replace=False))]
N = len(users)
print(f"S4 control. users with all four (window x half) groups: {N:,}")
print(f"targets from O3: S2-NN top-1 {TARGET_TOP1:.4f}, transported ECE {S2_ECE_TR:.4f}\n")

def xy(k): a, b = k.split("_"); return int(a), int(b)
def anchor(c): return xy(max(c.items(), key=lambda kv: (kv[1], kv[0]))[0])
def offs(c, a):
    tot = sum(c.values()); out = {}
    for k, v in c.items():
        x, y = xy(k); out[(x-a[0], y-a[1])] = v/tot
    return out
def jsd(p, q):
    keys = sorted(set(p) | set(q))
    return M.js(np.array([p.get(k,0.0) for k in keys]), np.array([q.get(k,0.0) for k in keys]))

P = {}
for u in users:
    e = {}
    for half in "MD":
        aR = anchor(U[u]["R"+half]["c"]); aT = anchor(U[u]["T"+half]["c"])
        e[half] = dict(aR=aR, aT=aT, qR=offs(U[u]["R"+half]["c"], aR),
                       qT=offs(U[u]["T"+half]["c"], aT))
    P[u] = e

d_own = np.array([jsd(P[u]["D"]["qR"], P[u]["D"]["qT"]) for u in users])
covD  = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)

# ---------------------------------------------------------------- the perturbation
DIRS = {u: rng.dirichlet(np.ones(len(P[u]["M"]["qT"]))) for u in users}
def perturb(u, lam):
    q = P[u]["M"]["qT"]; ks = list(q); w = DIRS[u]
    mix = {k: (1-lam)*q[k] + lam*float(w[i]) for i, k in enumerate(ks)}
    s = sum(mix.values())
    return {k: v/s for k, v in mix.items() if v > 1e-12}

def to_cells(q, a, total):
    c = {}
    for (dx, dy), m in q.items():
        x, y = a[0]+dx, a[1]+dy
        if 0 <= x < GRID and 0 <= y < GRID: c[f"{x}_{y}"] = c.get(f"{x}_{y}", 0.0)+m
    s = sum(c.values())
    return ({k: v*total/s for k, v in c.items()} if s > 0 else None)

def build(lams):
    """lams: array of per-user lambda. Returns the synthetic world."""
    Uc = {}
    for i, u in enumerate(users):
        q = perturb(u, float(lams[i]))
        c = to_cells(q, P[u]["M"]["aT"], U[u]["TM"]["n"])
        if c is None: continue
        Uc[u] = {"RM": U[u]["RM"], "RD": U[u]["RD"],
                 "TM": dict(n=U[u]["TM"]["n"], d=U[u]["TM"]["d"], n27=0, c=c, s=U[u]["TM"]["s"]),
                 "TD": U[u]["TD"]}          # <-- the D half is the REAL one, untouched
    return Uc

# ---------------------------------------------------------------- real world + calibrator
print("REAL world (no perturbation), to fit the calibrator that will be transported")
res0 = M.match(U, users, "RM", "TM")
cal = rng.permutation(N); ev = cal[N//2:]; cal = cal[:N//2]
k0 = M.calibrate(res0, covD, cal, ev)
REALCAL = dict(w=k0["w"], mu=k0["mu"], sd=k0["sd"])
print(f"  top-1 {res0['y'].mean():.4f}   ECE(eval) {M.ece(k0['p_platt'][ev], res0['y'][ev]):.4f}\n")

def run(Uc, label):
    sub = [u for u in users if u in Uc]
    ii = {u: j for j, u in enumerate(sub)}
    idx = np.array([i for i, u in enumerate(users) if u in ii])
    res = M.match(Uc, sub, "RM", "TM")
    cv = covD[idx]
    X = M.platt_features(res["topk"], cv)
    p_tr = M.logit_apply((X - REALCAL["mu"])/REALCAL["sd"], REALCAL["w"])
    evm = np.array([ii[users[i]] for i in ev if users[i] in ii])
    y = res["y"]
    # measured drift in this world, on the UNTOUCHED D half
    d_app = np.array([jsd(P[u]["D"]["qR"], P[u]["D"]["qT"]) for u in sub])
    return dict(label=label, sub=sub, idx=idx, y=y, p_tr=p_tr, ev=evm, d_app=d_app,
                top1=float(y.mean()), ece=float(M.ece(p_tr[evm], y[evm])))

# ---------------------------------------------------------------- S4-ACC: bisect global lam
print("S4-ACC: bisecting one global lambda so pool top-1 matches S2-NN")
lo, hi = 0.0, 1.0
for it in range(8):
    mid = 0.5*(lo+hi)
    r = run(build(np.full(N, mid)), f"probe lam={mid:.4f}")
    print(f"  lam {mid:.4f}  top-1 {r['top1']:.4f}")
    if r["top1"] > TARGET_TOP1: lo = mid
    else: hi = mid
LAM = 0.5*(lo+hi)
S4A = run(build(np.full(N, LAM)), "S4-ACC")
print(f"  chosen lambda {LAM:.4f}\n")

# ---------------------------------------------------------------- lambda sweep
print("sweep: fixed lambdas, to show the whole disruption-versus-ECE curve")
SWEEP = []
for lam in (0.2, 0.4, 0.6, 0.8, 1.0):
    r = run(build(np.full(N, lam)), f"lam={lam:.1f}")
    SWEEP.append(r)
    print(f"  lam {lam:.1f}  top-1 {r['top1']:.4f}  ECE_transported {r['ece']:.4f}")
print()

# ---------------------------------------------------------------- verdict
print("="*94)
print("PRIMARY: transported-calibrator ECE, against S2-NN")
print("="*94)
print(f"  {'condition':<12} {'n':>7} {'top-1':>8} {'ECE_transported':>17}")
print(f"  {'REAL':<12} {N:>7,} {res0['y'].mean():>8.4f} {M.ece(k0['p_platt'][ev],res0['y'][ev]):>17.4f}")
for r in (S4A,):
    print(f"  {r['label']:<12} {len(r['sub']):>7,} {r['top1']:>8.4f} {r['ece']:>17.4f}")
print(f"  {'S2-NN (O3)':<12} {'6,224':>7} {TARGET_TOP1:>8.4f} {S2_ECE_TR:>17.4f}")

e = S4A["ece"]
print(f"\n  ratio ECE(S4-ACC)/ECE(S2-NN) = {e/S2_ECE_TR:.3f}")
if abs(S4A["top1"] - TARGET_TOP1) > 0.01:
    print(f"  !! top-1 not matched within 0.01 ({S4A['top1']:.4f} vs {TARGET_TOP1:.4f}); "
          "the comparison is NOT at matched difficulty and the verdict is void")
elif e >= 0.8*S2_ECE_TR: print("  VERDICT: WITHDRAW -- generic disruption reproduces the failure")
elif e <= 0.5*S2_ECE_TR: print("  VERDICT: SURVIVES -- drift structure is doing the work")
else:                    print("  VERDICT: INCONCLUSIVE")

print("\nINTEGRITY: measured drift must be unchanged by construction")
for r in (S4A,):
    dd = np.abs(r["d_app"] - d_own[r["idx"]]).max()
    print(f"  {r['label']:<12} max|d_app - d_own| = {dd:.3e}  "
          f"{'OK' if dd < 1e-12 else '<<< BROKEN'}")

print("\nSECONDARY: residual on d_app within S4 (must be ~0; d_app is constant by construction)")
def z(v): return (v - v.mean())/max(v.std(), 1e-12)
for r in (S4A,):
    evm = r["ev"]; resid = r["y"][evm] - r["p_tr"][evm]
    X = np.column_stack([np.ones(len(evm)), z(r["d_app"][evm])])
    b = np.linalg.lstsq(X, resid, rcond=None)[0]
    print(f"  {r['label']:<12} b(d_app) = {b[1]:+.4f}")

print("\nper-bin gap by applied-drift quintile (S4-ACC), for comparison with O3 §7")
evm = S4A["ev"]; dv = S4A["d_app"][evm]
q = np.quantile(dv, [0, .2, .4, .6, .8, 1.0])
for k in range(5):
    m = (dv >= q[k]) & (dv <= q[k+1] if k == 4 else dv < q[k+1])
    print(f"  Q{k+1} n={m.sum():>5}  mean d {dv[m].mean():.3f}  "
          f"gap {S4A['y'][evm][m].mean()-S4A['p_tr'][evm][m].mean():+.3f}")

np.savez("results/s4_control.npz",
         **{f"S4ACC_{f}": S4A[f] for f in ("y", "p_tr", "d_app", "idx", "ev")},
         sweep_lam=np.array([0.2,0.4,0.6,0.8,1.0]),
         sweep_top1=np.array([r["top1"] for r in SWEEP]),
         sweep_ece=np.array([r["ece"] for r in SWEEP]),
         lam_acc=np.array([LAM]), d_own=d_own)
print("\nsaved results/s4_control.npz")
