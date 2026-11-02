"""Discharge premises P3 (drift varies at fixed gap), P6 (pool size) and P7 (observability).

Usage: sweeps.py <windows.tsv> <label> [pool_size] [seed]
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M

PROF, LABEL = sys.argv[1], sys.argv[2]
POOL = int(sys.argv[3]) if len(sys.argv) > 3 else 0
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261102
GAPS = [(0, "W3"), (15, "W2"), (30, "W1"), (45, "W0")]
TGT = "W4"

U = M.parse_windows(PROF)
need = {f"W{w}{h}" for w in range(5) for h in "MD"}
complete = sorted([u for u, g in U.items() if need <= set(g)], key=int)
rng = np.random.default_rng(SEED)
if POOL and POOL < len(complete):
    users = [complete[i] for i in sorted(rng.choice(len(complete), POOL, replace=False))]
else:
    users = complete
N = len(users)
print(f"=== {LABEL} ===")
print(f"users present in all 10 (window x half) groups: {len(complete):,}; pool used: {N:,}\n")

cov = {g: np.array([U[u][g]["n"] for u in users], float) for g in need}
perm = rng.permutation(N); cal, ev = perm[:N//2], perm[N//2:]

# ================================================================= P3
print("="*100)
print("P3 - does drift vary at FIXED gap, and is drift collinear with gap?")
print("="*100)
store = {}
rows = []
for g, rw in GAPS:
    res = M.match(U, users, rw + "M", TGT + "M")
    d = M.drift(U, users, rw + "D", TGT + "D")
    c = cov[rw + "D"] + cov[TGT + "D"]
    k = M.calibrate(res, c, cal, ev)
    store[g] = dict(res=res, d=d, cov=c, k=k)
    js = d["js_cells"]
    rows.append((g, rw, res["y"].mean(), js.mean(), js.std(), np.quantile(js, .1),
                 np.quantile(js, .9), M.spear(js, k["p_platt"]), M.ece(k["p_platt"][ev], res["y"][ev])))
print(f"  {'gap':>4} {'ref':>4} {'top-1':>7} {'mean JS':>8} {'sd JS':>7} {'p10':>6} {'p90':>6} "
      f"{'r(JS,p)':>8} {'ECE(ev)':>8}")
for r in rows:
    print(f"  {r[0]:>4} {r[1]:>4} {r[2]:>7.4f} {r[3]:>8.3f} {r[4]:>7.3f} {r[5]:>6.3f} {r[6]:>6.3f} "
          f"{r[7]:>+8.3f} {r[8]:>8.4f}")

allg = np.concatenate([[g]*N for g, _ in GAPS])
alljs = np.concatenate([store[g]["d"]["js_cells"] for g, _ in GAPS])
print(f"\n  pooled over all four gaps: r(JS drift, gap) = {M.spear(alljs, allg):+.3f}")
wi = np.mean([store[g]['d']['js_cells'].std() for g, _ in GAPS])
print(f"  mean WITHIN-gap sd of JS drift = {wi:.3f}; BETWEEN-gap sd of the means = "
      f"{np.std([r[3] for r in rows]):.3f}  -> within/between = {wi/max(np.std([r[3] for r in rows]),1e-9):.1f}x")
print(f"  at the primary g = 0 the gap is constant by construction, and JS drift still spans "
      f"[{store[0]['d']['js_cells'].min():.3f}, {store[0]['d']['js_cells'].max():.3f}]")
print("\n  day 27 exposure (Data Descriptor advises excluding it; it falls in W1, used at g = 30):")
for g, rw in GAPS:
    n27 = sum(U[u][rw+"M"]["n27"] + U[u][rw+"D"]["n27"] for u in users)
    tot = sum(U[u][rw+"M"]["n"] + U[u][rw+"D"]["n"] for u in users)
    print(f"    g={g:<3} ref={rw}: {n27:,} of {tot:,} records ({100*n27/tot:.2f}%) fall on day 27")

print("\n  nested fit WITHIN each gap stratum (does drift add beyond confidence + coverage?):")
print(f"  {'gap':>4} {'AUC c+cov':>10} {'AUC +drift':>11} {'dAUC':>7} {'LL c+cov':>9} {'LL +drift':>10}")
def auc(p, y):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1*(n1+1)/2)/(n1*n0))
def z(v): return (v - v.mean())/max(v.std(), 1e-12)
for g, rw in GAPS:
    s = store[g]; y = s["res"]["y"]
    base = np.column_stack([z(s["k"]["p_platt"]), z(np.log(s["cov"]))])
    full = np.column_stack([base, z(s["d"]["js_cells"])])
    o = []
    for X in (base, full):
        w = M.logit_fit(X[cal], y[cal]); p = M.logit_apply(X, w)
        o.append((auc(p[ev], y[ev]), -M.nll(p[ev], y[ev])))
    print(f"  {g:>4} {o[0][0]:>10.4f} {o[1][0]:>11.4f} {o[1][0]-o[0][0]:>+7.4f} "
          f"{o[0][1]:>9.4f} {o[1][1]:>10.4f}")

# ================================================================= P6
print("\n" + "="*100)
print("P6 - sensitivity to pool size (confidence is pool-dependent by construction)")
print("="*100)
print(f"  {'pool':>7} {'top-1':>7} {'ECE':>8} {'Brier':>8} {'T':>6} {'r(JS,p)':>8} "
      f"{'gap in top JS bin':>18} {'n in bin':>9}")
for ps in [1000, 2500, 5000, 10000, N]:
    if ps > N: continue
    idx = np.sort(rng.choice(N, ps, replace=False)) if ps < N else np.arange(N)
    su = [users[i] for i in idx]
    res = M.match(U, su, "W3M", TGT + "M")
    d = M.drift(U, su, "W3D", TGT + "D")
    c = np.array([U[u]["W3D"]["n"] + U[u][TGT+"D"]["n"] for u in su], float)
    p2 = rng.permutation(ps); ca, evp = p2[:ps//2], p2[ps//2:]
    k = M.calibrate(res, c, ca, evp)
    y, p = res["y"], k["p_platt"]
    m = d["js_cells"][evp] >= 0.85
    gp = (y[evp][m].mean() - p[evp][m].mean()) if m.sum() >= 30 else float("nan")
    print(f"  {ps:>7,} {y.mean():>7.4f} {M.ece(p[evp],y[evp]):>8.4f} {M.brier(p[evp],y[evp]):>8.4f} "
          f"{k['T']:>6.3f} {M.spear(d['js_cells'],p):>+8.3f} {gp:>18.3f} {m.sum():>9,}")

# ================================================================= P7
print("\n" + "="*100)
print("P7 - observability. Sensitivity to a per-user coverage filter at the primary g = 0.")
print("="*100)
s = store[0]; y = s["res"]["y"]; p = s["k"]["p_platt"]; js = s["d"]["js_cells"]; c = s["cov"]
cmin = np.minimum(cov["W3D"], cov[TGT+"D"])
print("  'coverage' = records in the DRIFT half; the filter uses the MINIMUM of the two windows,")
print("  so a user must be observed in both. Thresholds are reported as a curve, not chosen:")
print(f"  {'min recs >=':>12} {'users kept':>11} {'%':>6} {'top-1':>7} {'ECE':>8} "
      f"{'r(JS,p)':>8} {'r(JS,cov)':>10} {'top-bin gap':>12} {'n':>7}")
for th in [0, 20, 40, 60, 80, 120, 160]:
    m = cmin >= th
    if m.sum() < 200: continue
    me = m[ev]
    mb = me & (js[ev] >= 0.85)
    gp = (y[ev][mb].mean() - p[ev][mb].mean()) if mb.sum() >= 30 else float("nan")
    print(f"  {th:>12} {m.sum():>11,} {100*m.sum()/N:>6.1f} {y[m].mean():>7.4f} "
          f"{M.ece(p[ev][me],y[ev][me]):>8.4f} {M.spear(js[m],p[m]):>+8.3f} "
          f"{M.spear(js[m],c[m]):>+10.3f} {gp:>12.3f} {mb.sum():>7,}")
print("\n  entropy-change drift is the sample-size-sensitive component; its coverage")
print("  dependence at g = 0, before and after conditioning on coverage decile:")
de = s["d"]["d_entropy"]
dec = np.digitize(c, np.quantile(c, np.arange(.1, 1.0, .1)))
resid = de - np.array([de[dec == k].mean() if (dec == k).sum() else 0 for k in dec])
print(f"    r(d_entropy, coverage) raw = {M.spear(de, c):+.3f};  "
      f"within-decile residual vs coverage = {M.spear(resid, c):+.3f}")
print(f"    r(d_entropy, p) raw = {M.spear(de, p):+.3f};  residual = {M.spear(resid, p):+.3f}")

np.savez(f"results/sweeps_{LABEL}.npz",
         **{f"g{g}_y": store[g]["res"]["y"] for g, _ in GAPS},
         **{f"g{g}_p": store[g]["k"]["p_platt"] for g, _ in GAPS},
         **{f"g{g}_js": store[g]["d"]["js_cells"] for g, _ in GAPS},
         **{f"g{g}_cov": store[g]["cov"] for g, _ in GAPS}, cal=cal, ev=ev)
print(f"\nsaved results/sweeps_{LABEL}.npz")
