"""O3 stage 2: apply the pre-registered decision rule of B_SWAP_DESIGN.md §7."""
import pickle, numpy as np
from scipy import stats

R = pickle.load(open("results/o3_result.pkl", "rb"))
EDGES = [0.0, 0.45, 0.55, 0.65, 0.75, 0.85, 1.01]
rng = np.random.default_rng(20261102)

def ece(p, y, nb=15):
    e = 0.0
    for k in range(nb):
        m = (p >= k/nb) & (p < (k+1)/nb if k < nb-1 else p <= 1.0)
        if m.sum(): e += m.sum()/len(p)*abs(p[m].mean() - y[m].mean())
    return float(e)
def z(v): return (v - v.mean())/max(v.std(), 1e-12)
def ols(X, y):
    X = np.column_stack([np.ones(len(y)), X])
    return np.linalg.lstsq(X, y, rcond=None)[0]
def boot_ols(X, y, reps=2000):
    out = np.zeros((reps, X.shape[1]+1))
    for r in range(reps):
        i = rng.integers(0, len(y), len(y))
        out[r] = ols(X[i], y[i])
    return np.quantile(out, [.025, .975], axis=0)

print("="*104)
print("1. Gate and fidelity, on the pre-registered thresholds (admissible >= 50%, S2 delivery >= 60%)")
print("="*104)
print(f"  {'condition':<20} {'n':>7} {'admissible':>11} {'top-1':>7} {'ECE(ev)':>8} "
      f"{'mean d_app':>11} {'d_app/d_donor':>14} {'r(d_app,d_own)':>15}")
for k, d in R.items():
    da, do, dd, adm = d["da"], d["do"], d["dd"], d["adm"]
    deliv = (np.nanmean(da)/np.nanmean(dd)) if np.isfinite(dd).any() else float("nan")
    print(f"  {k:<20} {len(da):>7,} {100*adm.mean():>10.1f}% {d['y'].mean():>7.4f} "
          f"{ece(d['p'][d['ev']], d['y'][d['ev']]):>8.4f} {np.nanmean(da):>11.3f} "
          f"{deliv:>14.2f} {stats.spearmanr(da, do).statistic:>+15.3f}")

print("\n" + "="*104)
print("2. THE DECISION TEST. Calibration residual (y - p) on applied vs own drift, both")
print("   z-scored, OLS on ADMISSIBLE evaluation users. Bootstrap 95% CI, 2,000 resamples.")
print("="*104)
print(f"  {'condition':<20} {'n':>6} {'b(d_app)':>10} {'95% CI':>20} {'b(d_own)':>10} "
      f"{'95% CI':>20} {'VIF':>6}")
VERDICT = {}
for k, d in R.items():
    ev = d["ev"]; m = ev[d["adm"][ev]]
    if k.startswith(("REAL", "S3")):
        print(f"  {k:<20} {len(m):>6} -- not testable: d_app is identical to d_own by "
              f"construction (self/no swap)"); VERDICT[k] = "n/a"; continue
    if len(m) < 200:
        print(f"  {k:<20} {len(m):>6} -- too few admissible users to test"); VERDICT[k] = "no test"; continue
    da, do = d["da"][m], d["do"][m]
    if not (np.isfinite(da).all() and np.isfinite(do).all()):
        print(f"  {k:<20} {len(m):>6} -- non-finite drift"); VERDICT[k] = "no test"; continue
    r = abs(stats.spearmanr(da, do).statistic)
    vif = 1/max(1e-9, 1 - r**2)
    X = np.column_stack([z(da), z(do)])
    resid = d["y"][m] - d["p"][m]
    b = ols(X, resid); ci = boot_ols(X, resid)
    sig = lambda j: not (ci[0][j] <= 0 <= ci[1][j])
    print(f"  {k:<20} {len(m):>6} {b[1]:>+10.4f} [{ci[0][1]:>+8.4f},{ci[1][1]:>+8.4f}] "
          f"{b[2]:>+10.4f} [{ci[0][2]:>+8.4f},{ci[1][2]:>+8.4f}] {vif:>6.1f}")
    VERDICT[k] = ("app only" if sig(1) and not sig(2) else
                  "own only" if sig(2) and not sig(1) else
                  "both" if sig(1) and sig(2) else "neither")

print("\n" + "="*104)
print("3. Drift-stratified calibration gap on ABSOLUTE d_app bins, admissible evaluation users")
print("="*104)
for k, d in R.items():
    ev = d["ev"]; m = ev[d["adm"][ev]]
    if len(m) < 200: continue
    print(f"\n  {k}  (n = {len(m):,} admissible of {len(ev):,} evaluation users)")
    print(f"    {'d_app bin':<12} {'n':>6} {'acc':>7} {'mean p':>8} {'gap':>8} {'ECE':>7}")
    for j in range(len(EDGES)-1):
        lo, hi = EDGES[j], EDGES[j+1]
        b = m[(d["da"][m] >= lo) & (d["da"][m] < hi)]
        if len(b) < 30: continue
        y, p = d["y"][b], d["p"][b]
        print(f"    {lo:.2f}-{min(hi,1.0):.2f}    {len(b):>6} {y.mean():>7.3f} {p.mean():>8.3f} "
              f"{y.mean()-p.mean():>+8.3f} {ece(p,y):>7.4f}")

print("\n" + "="*104)
print("4. Verdict against the rule fixed in B_SWAP_DESIGN.md §7")
print("="*104)
for k in R:
    adm = 100*R[k]["adm"].mean()
    print(f"  {k:<20} admissible {adm:>5.1f}%   residual tracks: {VERDICT.get(k,'-')}")

# ============================================================ AMENDED (exploratory)
print("\n" + "="*104)
print("5. AMENDED ANALYSIS (exploratory, A_PREMISES s8 2026-11-02). Same decision test, but on the")
print("   probability of the REAL world's calibrator TRANSPORTED unchanged into the swapped world.")
print("   The pre-registered refit absorbed the applied drift; this is where O1's effect lives.")
print("="*104)
have = all("p_tr" in d for d in R.values())
if not have:
    print("  o3_result.pkl predates the amendment; rerun scripts/o3_swap.py")
else:
    print(f"  {'condition':<20} {'ECE refit':>10} {'ECE transp':>11} {'n':>6} "
          f"{'b(d_app)':>10} {'95% CI':>20} {'b(d_own)':>10} {'95% CI':>20}")
    for k, d in R.items():
        ev = d["ev"]; m = ev[d["adm"][ev]]
        e1 = ece(d["p"][ev], d["y"][ev]); e2 = ece(d["p_tr"][ev], d["y"][ev])
        if k.startswith(("REAL", "S3")) or len(m) < 200:
            print(f"  {k:<20} {e1:>10.4f} {e2:>11.4f} {len(m):>6}   (not testable)")
            continue
        da, do = d["da"][m], d["do"][m]
        X = np.column_stack([z(da), z(do)])
        resid = d["y"][m] - d["p_tr"][m]
        b = ols(X, resid); ci = boot_ols(X, resid)
        print(f"  {k:<20} {e1:>10.4f} {e2:>11.4f} {len(m):>6} {b[1]:>+10.4f} "
              f"[{ci[0][1]:>+8.4f},{ci[1][1]:>+8.4f}] {b[2]:>+10.4f} [{ci[0][2]:>+8.4f},{ci[1][2]:>+8.4f}]")

    print("\n  Drift-stratified gap under the TRANSPORTED calibrator, admissible evaluation users")
    for k, d in R.items():
        if k.startswith(("REAL", "S3")): continue
        ev = d["ev"]; m = ev[d["adm"][ev]]
        if len(m) < 200: continue
        print(f"\n  {k}  (n = {len(m):,})")
        print(f"    {'d_app bin':<12} {'n':>6} {'acc':>7} {'mean p':>8} {'gap':>8} {'ECE':>7}")
        for j in range(len(EDGES)-1):
            lo, hi = EDGES[j], EDGES[j+1]
            b = m[(d["da"][m] >= lo) & (d["da"][m] < hi)]
            if len(b) < 30: continue
            y, p = d["y"][b], d["p_tr"][b]
            print(f"    {lo:.2f}-{min(hi,1.0):.2f}    {len(b):>6} {y.mean():>7.3f} {p.mean():>8.3f} "
                  f"{y.mean()-p.mean():>+8.3f} {ece(p,y):>7.4f}")
