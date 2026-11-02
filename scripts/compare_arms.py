"""Cross-arm comparison of calibration, on an ABSOLUTE drift scale.

The per-arm composite in calibrate.py is rank-based, so its quintiles are not comparable
between arms. This script bins on visit-distribution JS computed on the held-out day half
-- a fixed [0,1] scale -- and compares identical absolute bins across arms.
"""
import numpy as np
A = np.load("results/cal_armA_dataset2.npz")
B = np.load("results/cal_armB_dataset1.npz")

def ece(p, y, nb=15):
    e = 0.0
    for k in range(nb):
        m = (p >= k/nb) & (p < (k+1)/nb if k < nb-1 else p <= 1.0)
        if m.sum(): e += m.sum()/len(p) * abs(p[m].mean() - y[m].mean())
    return float(e)

print("Absolute drift distribution (visit-distribution JS, held-out half), all users")
print(f"  {'arm':<6} {'n':>7} {'mean':>7} {'p10':>7} {'p25':>7} {'p50':>7} {'p75':>7} {'p90':>7} {'p99':>7}")
for nm, d in (("A", A), ("B", B)):
    v = d["dec_js_cells"]
    q = np.quantile(v, [.10, .25, .50, .75, .90, .99])
    print(f"  {nm:<6} {len(v):>7,} {v.mean():>7.3f}" + "".join(f" {x:>7.3f}" for x in q))

EDGES = [0.0, 0.45, 0.55, 0.65, 0.75, 0.85, 1.01]
print(f"\nCommon absolute bins on JS drift. 'gap' = accuracy - mean predicted probability.")
print("Evaluation users only; the calibrator was fitted per arm on that arm's calibration users.")
print(f"\n  {'bin':<12} {'armA n':>7} {'A acc':>7} {'A mean p':>9} {'A gap':>8} {'A ECE':>7} | "
      f"{'armB n':>7} {'B acc':>7} {'B mean p':>9} {'B gap':>8} {'B ECE':>7}")
rows = []
for k in range(len(EDGES) - 1):
    lo, hi = EDGES[k], EDGES[k+1]
    out = [f"  {lo:.2f}-{min(hi,1.0):.2f}  "]
    rec = dict(lo=lo, hi=hi)
    for nm, d in (("A", A), ("B", B)):
        ev = d["ev"]; v = d["dec_js_cells"][ev]; p = d["p_platt"][ev]; y = d["y"][ev]
        m = (v >= lo) & (v < hi)
        if m.sum() < 30:
            out.append(f" {m.sum():>7} {'-':>7} {'-':>9} {'-':>8} {'-':>7}")
            rec[nm] = None
        else:
            g = y[m].mean() - p[m].mean()
            out.append(f" {m.sum():>7} {y[m].mean():>7.3f} {p[m].mean():>9.3f} {g:>+8.3f} {ece(p[m],y[m]):>7.4f}")
            rec[nm] = (int(m.sum()), float(y[m].mean()), float(p[m].mean()), float(g), ece(p[m], y[m]))
        if nm == "A": out.append(" |")
    print("".join(out)); rows.append(rec)

print("\nSame bins, but each arm's calibrator REPLACED by the other arm's, to separate")
print("'the shock changed behaviour' from 'the calibrator was fitted on different data'.")
def apply_w(d, w, mu_sd):
    s1, s2, topk, cov = d["s1"], d["s2"], d["topk"], d["cov"]
    X = np.column_stack([s1, s1-s2, s2-topk[:,9], s1-topk[:,:64].mean(1), np.log(cov)])
    mu, sd = mu_sd
    Xz = (X-mu)/sd
    return 1/(1+np.exp(-np.clip(np.column_stack([np.ones(len(s1)), Xz])@w, -30, 30)))
def musd(d):
    cal = d["cal"]
    X = np.column_stack([d["s1"], d["s1"]-d["s2"], d["s2"]-d["topk"][:,9],
                         d["s1"]-d["topk"][:,:64].mean(1), np.log(d["cov"])])
    return X[cal].mean(0), X[cal].std(0)
pA_fromB = apply_w(A, B["w"], musd(B))
print(f"  {'bin':<12} {'n':>7} {'A acc':>7} {'p (B-cal)':>10} {'gap':>8} {'ECE':>7}")
for k in range(len(EDGES)-1):
    lo, hi = EDGES[k], EDGES[k+1]
    ev = A["ev"]; v = A["dec_js_cells"][ev]; m = (v>=lo)&(v<hi)
    if m.sum() < 30: continue
    p, y = pA_fromB[ev][m], A["y"][ev][m]
    print(f"  {lo:.2f}-{min(hi,1.0):.2f}   {m.sum():>7} {y.mean():>7.3f} {p.mean():>10.3f} "
          f"{y.mean()-p.mean():>+8.3f} {ece(p,y):>7.4f}")

print("\nOverall, evaluation users:")
for nm, d in (("A (shock)", A), ("B (no shock)", B)):
    ev = d["ev"]
    print(f"  {nm:<14} top-1 {d['y'][ev].mean():.4f}  ECE {ece(d['p_platt'][ev], d['y'][ev]):.4f}  "
          f"Brier {np.mean((d['p_platt'][ev]-d['y'][ev])**2):.4f}  T {float(d['T']):.3f}")

# ---------------------------------------------------------------- ECE decomposition
print("\n" + "="*92)
print("Aggregate ECE vs drift-conditional ECE: does aggregate calibration hide the failure?")
print("="*92)
def cond_ece(d, p=None, edges=EDGES):
    ev = d["ev"]; v = d["dec_js_cells"][ev]
    p = (d["p_platt"] if p is None else p)[ev]; y = d["y"][ev]
    tot, n = 0.0, 0
    for k in range(len(edges)-1):
        m = (v >= edges[k]) & (v < edges[k+1])
        if m.sum() >= 30: tot += m.sum()*ece(p[m], y[m]); n += m.sum()
    return tot/n
for nm, d in (("A (shock)", A), ("B (no shock)", B)):
    ev = d["ev"]
    agg = ece(d["p_platt"][ev], d["y"][ev]); con = cond_ece(d)
    print(f"  {nm:<14} aggregate ECE {agg:.4f}   drift-conditional ECE {con:.4f}   "
          f"understated by {con/agg:.1f}x")
agg = ece(pA_fromB[A["ev"]], A["y"][A["ev"]])
print(f"  {'A, B-calibrator':<14} aggregate ECE {agg:.4f}   drift-conditional ECE "
      f"{cond_ece(A, pA_fromB):.4f}   understated by {cond_ece(A,pA_fromB)/agg:.1f}x")

# ---------------------------------------------------------------- bootstrap
print("\n" + "="*92)
print("Bootstrap 95% CI (2,000 resamples of evaluation users) for the top drift bin, JS >= 0.85")
print("="*92)
rng = np.random.default_rng(7)
def boot(v, y, p, reps=2000):
    g = np.zeros(reps); e = np.zeros(reps)
    for r in range(reps):
        i = rng.integers(0, len(y), len(y))
        g[r] = y[i].mean() - p[i].mean(); e[r] = ece(p[i], y[i])
    return np.quantile(g, [.025, .975]), np.quantile(e, [.025, .975])
for nm, d, p in (("A, own calibrator", A, None), ("A, B calibrator", A, pA_fromB), ("B, own calibrator", B, None)):
    ev = d["ev"]; v = d["dec_js_cells"][ev]; m = v >= 0.85
    pp = (d["p_platt"] if p is None else p)[ev][m]; yy = d["y"][ev][m]
    (gl, gh), (el, eh) = boot(v[m], yy, pp)
    print(f"  {nm:<20} n={m.sum():>5}  gap {yy.mean()-pp.mean():+.3f} [{gl:+.3f}, {gh:+.3f}]   "
          f"ECE {ece(pp,yy):.3f} [{el:.3f}, {eh:.3f}]")

print("\nDrift-score association in both arms (Spearman, all users):")
def spear(a, b):
    """Tie-corrected Spearman (midranks). An argsort-of-argsort implementation is WRONG
    whenever either variable has ties -- the outcome y is binary, so it is all ties, and the
    error attenuated every r(., y) in the 2026-11-02 runs by 0.10 to 0.19."""
    from scipy import stats
    return float(stats.spearmanr(a, b).statistic)
for nm, d in (("A (shock)", A), ("B (no shock)", B)):
    v = d["dec_js_cells"]
    print(f"  {nm:<14} r(JS drift, p_platt) = {spear(v, d['p_platt']):+.3f}   "
          f"r(JS drift, y) = {spear(v, d['y']):+.3f}   mean JS = {v.mean():.3f}")
