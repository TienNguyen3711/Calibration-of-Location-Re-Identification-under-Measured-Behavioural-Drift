"""The atypicality baseline (R1_ATYPICALITY_RULING.md section 6).

WHY THIS EXISTS. Yuksekgonul, Zhang, Zou & Guestrin, "Beyond Confidence: Reliable Models
Should Also Consider Atypicality", NeurIPS 2023 (DOI 10.52202/075280-1669), already report
expected calibration error BY PER-SAMPLE ATYPICALITY QUANTILE (their Table 1), where

    a_X(x) = - max_y log P(X = x | Y = y),   "with respect to the training distribution".

So the SHAPE of our O1 section 4b analysis is theirs. What is supposed to be ours is that the
stratifying covariate is a WITHIN-SUBJECT CHANGE between two periods, computed on records the
matcher never consumes. A reviewer who knows that paper asks the obvious question:

    compute the plain rarity of the TARGET window alone -- one window, no reference window,
    no parity split -- and stratify by that instead. Does the change measure explain anything
    that rarity does not?

If it does not, delta buys nothing and the paper has no contribution.

PRE-REGISTERED, written before the run (ruling section 6):
  (a) both predict and delta retains an independent coefficient -> contribution stands as stated;
  (b) atypicality absorbs delta -> the paper shrinks to "Yuksekgonul transfers to privacy
      attacks" and must be rewritten as that;
  (c) delta absorbs atypicality -> strongest outcome, to be stated most cautiously.

PRIMARY ENDPOINT, fixed here before seeing any number: dAUC(delta | base + atypicality) on
EVALUATION users, where base = [z(p_platt), z(log coverage)] -- identical to the C2b nested fit
in sweeps.py, so the number is directly comparable to the published +0.026...+0.030.
DECISION RULE:
  (b) if its bootstrap 95% CI contains 0 in Arm A;
  (c) if it excludes 0 AND dAUC(atypicality | base + delta) contains 0;
  (a) otherwise, provided dAUC(atypicality | base) also excludes 0.
SECONDARY, the arm contrast: delta's increment is ~0 in Arm B and +0.026...+0.030 in Arm A. If
atypicality shows the SAME arm contrast then it is shock-sensitive too, and the argument that
"a static rarity measure cannot know a shock happened" is WRONG and must be withdrawn.

THREE ESTIMATORS, spanning their menu. All are computed on the TARGET window, MATCHER half
(TM) -- i.e. on the records the matcher does consume, which is the whole point of the baseline --
with the "training distribution" being the population of REFERENCE-window matcher halves (RM).

  a_pop   model-free multinomial surprisal: mean negative log population probability of this
          user's visited cells. Needs no matcher, no embedding, no scores. This is their
          Definition 2.1 collapsed to a single class, -log P(x).
  a_knn   their KNN estimator: mean distance to the k = 5 nearest reference profiles in the
          predictor's representation, = 1 - mean(top-5 cosine). Read off the stored top-64
          scores, so it is EXACTLY the matcher's own neighbourhood.
  a_maha  their vision estimator: a single Gaussian with tied covariance fitted on the
          predictor's representation (truncated SVD to d dims), Mahalanobis^2 of the target
          vector from the reference-population mean. One "class", there being no class label.

Nothing is refitted: the matcher, the calibrator and delta are read from results/cal_<label>.npz,
the record of the 2026-11-02 Arm A/B run. User order is reconstructed from the profile file by
calibrate.py's own rule and then CHECKED against the stored coverage and js_cells arrays; the
script aborts if they do not match to within floating-point tolerance.

Arm B was run at Arm A's pool size (a random subsample of 24,517 from 99,549, seed 20261102), so
`pool` must repeat calibrate.py's own draw -- the first draw from a fresh default_rng(seed) --
or the user order will not line up. The alignment check below is what catches a wrong draw.

Usage: atyp_baseline.py <profiles.tsv> <results/cal_LABEL.npz> <label> [pool] [svd_dim] [seed]
"""
import sys, math, numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
sys.path.insert(0, "scripts")
import matchlib as M

PROF, NPZ, LABEL = sys.argv[1], sys.argv[2], sys.argv[3]
POOL = int(sys.argv[4]) if len(sys.argv) > 4 else 0
SVD_D = int(sys.argv[5]) if len(sys.argv) > 5 else 64
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 20261102
NB = 2000  # bootstrap resamples

print(f"=== atypicality baseline: {LABEL} ===")
print(f"profiles {PROF}\nrecord   {NPZ}\npool     {POOL or 'all'}  "
      f"svd_dim {SVD_D}  seed {SEED}  bootstrap {NB}\n")

# ---------------------------------------------------------------- load, align, verify
U = {}
for line in open(PROF):
    f = line.rstrip("\n").split("\t")
    cells = {}
    for tok in f[5][1:].split():
        k, v = tok.rsplit(":", 1); cells[k] = int(v)
    slots = np.zeros(48)
    for tok in f[6][1:].split():
        k, v = tok.rsplit(":", 1); slots[int(k)] = int(v)
    U.setdefault(f[1], {})[f[2]] = dict(n=int(f[3]), d=int(f[4]), c=cells, s=slots)

complete = sorted([u for u, g in U.items() if {"RM", "RD", "TM", "TD"} <= set(g)], key=int)
if POOL and POOL < len(complete):
    # repeat calibrate.py's draw EXACTLY: the first draw from a fresh default_rng(SEED)
    pick = np.random.default_rng(SEED).choice(len(complete), POOL, replace=False)
    users = [complete[i] for i in sorted(pick)]
    print(f"pool: random subsample of {POOL:,} from {len(complete):,} complete users (seed {SEED})")
else:
    users = complete
    print(f"pool: all {len(complete):,} complete users")
N = len(users)
Z = np.load(NPZ)
y, p_platt, cov, comp = Z["y"], Z["p_platt"], Z["cov"], Z["comp"]
topk, cal, ev = Z["topk"], Z["cal"], Z["ev"]
js_cells = Z["dec_js_cells"]
assert len(y) == N, f"pool mismatch: profiles give {N}, record has {len(y)}"

cov_re = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)
js_re = np.array([M.drift_one(U[u]["RD"], U[u]["TD"])["js_cells"] for u in users])
dc, dj = np.abs(cov_re - cov).max(), np.abs(js_re - js_cells).max()
print(f"alignment check  max|d coverage| = {dc:.3e}   max|d js_cells| = {dj:.3e}")
assert dc == 0 and dj < 1e-12, "user order does not reproduce the stored run -- ABORT"
print(f"pool {N:,} users   cal {len(cal):,}   eval {len(ev):,}   top-1 {y.mean():.4f}\n")

# ---------------------------------------------------------------- a_pop (model-free)
pop = {}
for u in users:
    for k, v in U[u]["RM"]["c"].items(): pop[k] = pop.get(k, 0) + v
Vp = len(pop); tot = sum(pop.values())
# add-one over the observed reference vocabulary plus one bucket for anything unseen
den = tot + Vp + 1
logq = {k: math.log((v + 1)/den) for k, v in pop.items()}
log_unseen = math.log(1/den)
a_pop = np.zeros(N)
for i, u in enumerate(users):
    c = U[u]["TM"]["c"]; n = sum(c.values())
    a_pop[i] = -sum(v*logq.get(k, log_unseen) for k, v in c.items())/n
print(f"a_pop   reference vocabulary {Vp:,} cells, {tot:,} records; "
      f"mean {a_pop.mean():.3f} sd {a_pop.std():.3f} range [{a_pop.min():.2f}, {a_pop.max():.2f}]")

# ---------------------------------------------------------------- a_knn (their KNN estimator)
KNN = 5
a_knn = 1.0 - topk[:, :KNN].mean(1)
print(f"a_knn   1 - mean(top-{KNN} cosine); mean {a_knn.mean():.4f} sd {a_knn.std():.4f}")

# ---------------------------------------------------------------- a_maha (their vision estimator)
# rebuild the predictor's representation exactly as matchlib.match does
vocab = {}
for u in users:
    for g in ("RM", "TM"):
        for k in U[u][g]["c"]: vocab.setdefault(k, len(vocab))
V = len(vocab)
def build(grp):
    rows, cols, vals = [], [], []
    for i, u in enumerate(users):
        c = U[u][grp]["c"]; t = sum(c.values())
        for k, v in c.items():
            rows.append(i); cols.append(vocab[k]); vals.append(math.sqrt(v/t))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(N, V))
R, T = build("RM"), build("TM")
df = np.asarray((R > 0).sum(0)).ravel() + np.asarray((T > 0).sum(0)).ravel()
Dg = sparse.diags(np.log((2*N + 1)/(df + 1)))
R, T = R @ Dg, T @ Dg
def l2(A):
    n = np.sqrt(np.asarray(A.multiply(A).sum(1)).ravel()); n[n == 0] = 1
    return sparse.diags(1/n) @ A
R, T = l2(R), l2(T)
rng = np.random.default_rng(SEED)
Ur, Sr, Vt = svds(R.astype(np.float64), k=SVD_D, random_state=SEED)
o = np.argsort(-Sr); Vt = Vt[o]
PR, PT = np.asarray(R @ Vt.T), np.asarray(T @ Vt.T)
mu = PR.mean(0)
C = np.cov(PR - mu, rowvar=False) + 1e-9*np.eye(SVD_D)
Ci = np.linalg.inv(C)
dv = PT - mu
a_maha = np.einsum("ij,jk,ik->i", dv, Ci, dv)
print(f"a_maha  single Gaussian, tied cov, d = {SVD_D} (SVD captures "
      f"{(Sr[o]**2).sum()/ (R.multiply(R)).sum():.3f} of reference energy); "
      f"mean {a_maha.mean():.1f} median {np.median(a_maha):.1f}\n")

ATYP = {"a_pop": a_pop, "a_knn": a_knn, "a_maha": a_maha}

# ---------------------------------------------------------------- helpers
def auc(p, yv):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = yv.sum(); n0 = len(yv) - n1
    return float((r[yv == 1].sum() - n1*(n1 + 1)/2)/(n1*n0))
def z(v): return (v - v.mean())/max(v.std(), 1e-12)
def fitpred(X):
    w = M.logit_fit(X[cal], y[cal]); return M.logit_apply(X, w)

print("="*100)
print("correlations (tie-corrected Spearman, all users)")
print("="*100)
print(f"  {'':<8} {'r(.,y)':>9} {'r(.,p)':>9} {'r(.,cov)':>9} {'r(.,js)':>9} {'r(.,comp)':>10}")
for k, a in ATYP.items():
    print(f"  {k:<8} {M.spear(a,y):>+9.3f} {M.spear(a,p_platt):>+9.3f} {M.spear(a,cov):>+9.3f} "
          f"{M.spear(a,js_cells):>+9.3f} {M.spear(a,comp):>+10.3f}")
print(f"  {'js_cells':<8} {M.spear(js_cells,y):>+9.3f} {M.spear(js_cells,p_platt):>+9.3f} "
      f"{M.spear(js_cells,cov):>+9.3f} {'-':>9} {M.spear(js_cells,comp):>+10.3f}")
print("\n  pairwise among the atypicality estimators:")
ks = list(ATYP)
for i in range(len(ks)):
    for j in range(i+1, len(ks)):
        print(f"    r({ks[i]}, {ks[j]}) = {M.spear(ATYP[ks[i]], ATYP[ks[j]]):+.3f}")

# ---------------------------------------------------------------- T1: their Table 1 shape
print("\n" + "="*100)
print("T1 - reliability by ATYPICALITY quintile (the shape of their Table 1), eval users")
print("     then the same table by DRIFT quintile, for comparison (this is our O1 4b)")
print("="*100)
def strat(v, name):
    q = np.quantile(v[ev], [0, .2, .4, .6, .8, 1.0])
    print(f"\n  by {name}")
    print(f"  {'quintile':<9} {'n':>6} {'mean val':>10} {'accuracy':>9} {'mean p':>8} "
          f"{'gap':>8} {'ECE':>7} {'Brier':>7}")
    for k in range(5):
        m = (v[ev] >= q[k]) & (v[ev] <= q[k+1] if k == 4 else v[ev] < q[k+1])
        pe, ye = p_platt[ev][m], y[ev][m]
        print(f"  Q{k+1:<8} {m.sum():>6} {v[ev][m].mean():>10.3f} {ye.mean():>9.3f} "
              f"{pe.mean():>8.3f} {ye.mean()-pe.mean():>+8.3f} {M.ece(pe,ye):>7.4f} "
              f"{M.brier(pe,ye):>7.4f}")
for k, a in ATYP.items(): strat(a, k)
strat(js_cells, "js_cells (drift, ours)")
strat(comp, "comp (composite drift, ours)")

# ---------------------------------------------------------------- T2: the primary endpoint
print("\n" + "="*100)
print("T2 - PRIMARY. Nested dAUC on eval users. base = [z(p_platt), z(log coverage)],")
print("     identical to the C2b fit in sweeps.py. Bootstrap: eval users resampled, the")
print("     coefficients stay fitted on the (unresampled) calibration half.")
print("="*100)
base = np.column_stack([z(p_platt), z(np.log(cov))])
A3 = np.column_stack([z(a_pop), z(a_knn), z(a_maha)])
D1 = z(js_cells)[:, None]
MODELS = {
    "base (conf + coverage)":        base,
    "base + a_pop":                  np.column_stack([base, z(a_pop)]),
    "base + a_knn":                  np.column_stack([base, z(a_knn)]),
    "base + a_maha":                 np.column_stack([base, z(a_maha)]),
    "base + atyp (all 3)":           np.column_stack([base, A3]),
    "base + drift (js_cells)":       np.column_stack([base, D1]),
    "base + atyp + drift":           np.column_stack([base, A3, D1]),
}
P = {k: fitpred(X) for k, X in MODELS.items()}
boot = rng.integers(0, len(ev), size=(NB, len(ev)))
AU = {}
for k, p in P.items():
    pe, ye = p[ev], y[ev]
    bs = np.array([auc(pe[b], ye[b]) for b in boot])
    AU[k] = (auc(pe, ye), bs)
print(f"\n  {'model':<26} {'AUC':>8} {'NLL':>8}")
for k in MODELS:
    print(f"  {k:<26} {AU[k][0]:>8.4f} {M.nll(P[k][ev], y[ev]):>8.4f}")

def dauc(k1, k0):
    d = AU[k1][0] - AU[k0][0]
    bd = AU[k1][1] - AU[k0][1]
    lo, hi = np.quantile(bd, [0.025, 0.975])
    return d, lo, hi
print(f"\n  {'increment':<46} {'dAUC':>8}  {'95% CI':>18}")
TESTS = [
    ("atypicality (all 3) | base",            "base + atyp (all 3)",  "base (conf + coverage)"),
    ("drift | base   [the published C2b number]", "base + drift (js_cells)", "base (conf + coverage)"),
    ("*** drift | base + atypicality  PRIMARY", "base + atyp + drift",  "base + atyp (all 3)"),
    ("atypicality | base + drift",             "base + atyp + drift",  "base + drift (js_cells)"),
    ("a_pop | base",                           "base + a_pop",         "base (conf + coverage)"),
    ("a_knn | base",                           "base + a_knn",         "base (conf + coverage)"),
    ("a_maha | base",                          "base + a_maha",        "base (conf + coverage)"),
]
res = {}
for nm, k1, k0 in TESTS:
    d, lo, hi = dauc(k1, k0); res[nm] = (d, lo, hi)
    star = "" if lo <= 0 <= hi else "  CI excludes 0"
    print(f"  {nm:<46} {d:>+8.4f}  [{lo:>+7.4f},{hi:>+7.4f}]{star}")

# ---------------------------------------------------------------- T3: calibration residual
print("\n" + "="*100)
print("T3 - the calibration residual (y - p_platt) on eval users, OLS on z-scored covariates.")
print("     A negative coefficient means: more of this -> the calibrator is overconfident.")
print("="*100)
r = y[ev] - p_platt[ev]
def ols(cols, names):
    X = np.column_stack([np.ones(len(r))] + [z(c[ev]) for c in cols])
    b = np.linalg.lstsq(X, r, rcond=None)[0]
    bb = np.array([np.linalg.lstsq(X[bi], r[bi], rcond=None)[0] for bi in boot])
    print(f"\n  model: {' + '.join(names)}")
    for i, nm in enumerate(["(intercept)"] + names):
        lo, hi = np.quantile(bb[:, i], [0.025, 0.975])
        star = "" if lo <= 0 <= hi else "  *"
        print(f"    {nm:<14} {b[i]:>+8.4f}  [{lo:>+7.4f},{hi:>+7.4f}]{star}")
ols([a_pop], ["a_pop"]); ols([a_knn], ["a_knn"]); ols([a_maha], ["a_maha"])
ols([js_cells], ["js_cells"])
ols([a_pop, a_knn, a_maha], ["a_pop", "a_knn", "a_maha"])
ols([a_pop, a_knn, a_maha, js_cells], ["a_pop", "a_knn", "a_maha", "js_cells"])

np.savez(f"results/atyp_{LABEL}.npz", a_pop=a_pop, a_knn=a_knn, a_maha=a_maha,
         **{"dauc_" + k.replace(" ", "_"): np.array(v) for k, v in res.items()})
print(f"\nsaved results/atyp_{LABEL}.npz")
