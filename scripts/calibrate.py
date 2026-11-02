"""Give the matcher M a probability output, then measure its calibration (O1).

Two confidence outputs are produced and compared:

  softmax(T)  P(top-1 correct) = exp(s_1/T) / sum_v exp(s_v/T), the EXACT closed-world
              posterior over the pool under a Gibbs model of the score. One free
              parameter T, chosen on a grid by NLL on CALIBRATION users only.
  platt       logistic regression on score order statistics, fitted on CALIBRATION
              users only. No drift feature enters either calibrator, by design: O1 asks
              whether a drift-blind calibrator stays calibrated across drift strata.

Users are split 50/50 into calibration and evaluation sets. The candidate POOL is the
full pool in both cases -- the split decides whose metrics are reported, not who can be
matched. Pool size is an argument so arms can be compared at equal pool size (P6).

Drift uses the odd-day half only (P2, see results/P2_decoupling.md).

Usage: calibrate.py <profiles.tsv> <label> [pool_size] [seed]
"""
import sys, math, json, numpy as np
from scipy import sparse

PROF = sys.argv[1]
LABEL = sys.argv[2]
POOL = int(sys.argv[3]) if len(sys.argv) > 3 else 0
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261102
TGRID = [0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3]

def parse(path):
    U = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        cells = {}
        for tok in f[5][1:].split():
            k, v = tok.rsplit(":", 1); cells[k] = int(v)
        slots = np.zeros(48)
        for tok in f[6][1:].split():
            k, v = tok.rsplit(":", 1); slots[int(k)] = int(v)
        U.setdefault(f[1], {})[f[2]] = dict(n=int(f[3]), d=int(f[4]), c=cells, s=slots)
    return U

U = parse(PROF)
complete = sorted([u for u, g in U.items() if {"RM","RD","TM","TD"} <= set(g)], key=int)
rng = np.random.default_rng(SEED)
if POOL and POOL < len(complete):
    users = [complete[i] for i in sorted(rng.choice(len(complete), POOL, replace=False))]
    note = f"random subsample of {POOL:,} from {len(complete):,} complete users (seed {SEED})"
else:
    users = complete
    note = f"all {len(complete):,} complete users"
N = len(users)
print(f"=== {LABEL} ===")
print(f"pool: {note}\n")

# ------------------------------------------------------------------ matcher, M half only
vocab = {}
for u in users:
    for g in ("RM", "TM"):
        for k in U[u][g]["c"]: vocab.setdefault(k, len(vocab))
V = len(vocab)

def build(grp):
    rows, cols, vals = [], [], []
    for i, u in enumerate(users):
        c = U[u][grp]["c"]; tot = sum(c.values())
        for k, v in c.items():
            rows.append(i); cols.append(vocab[k]); vals.append(math.sqrt(v / tot))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(N, V))

R, T = build("RM"), build("TM")
df = np.asarray((R > 0).sum(0)).ravel() + np.asarray((T > 0).sum(0)).ravel()
Dg = sparse.diags(np.log((2 * N + 1) / (df + 1)))
R, T = R @ Dg, T @ Dg
def l2(M):
    n = np.sqrt(np.asarray(M.multiply(M).sum(1)).ravel()); n[n == 0] = 1
    return sparse.diags(1 / n) @ M
R, T = l2(R), l2(T)
Rt = R.T.tocsc()

K = 64
topk = np.zeros((N, K)); strue = np.zeros(N); rank = np.zeros(N, int); top = np.zeros(N, int)
lse = np.zeros((N, len(TGRID)))
B = 512
for a in range(0, N, B):
    b = min(a + B, N)
    S = np.asarray((T[a:b] @ Rt).todense())
    idx = np.argpartition(-S, K, axis=1)[:, :K]
    vals = np.take_along_axis(S, idx, 1)
    order = np.argsort(-vals, axis=1)
    topk[a:b] = np.take_along_axis(vals, order, 1)
    top[a:b] = np.take_along_axis(idx, order, 1)[:, 0]
    ar = np.arange(a, b)
    strue[a:b] = S[np.arange(b - a), ar]
    rank[a:b] = (S > strue[a:b, None]).sum(1) + 1
    mx = S.max(1, keepdims=True)
    for ti, t in enumerate(TGRID):
        lse[a:b, ti] = (mx[:, 0] / t) + np.log(np.exp(S / t - mx / t).sum(1))
y = (top == np.arange(N)).astype(float)
s1, s2 = topk[:, 0], topk[:, 1]
print(f"matcher: top-1 {y.mean():.4f}  top-10 {(rank<=10).mean():.4f}  median rank {np.median(rank):.0f}")

# ------------------------------------------------------------------ drift, D half only
def js(p, q):
    p, q = p / p.sum(), q / q.sum(); m = 0.5 * (p + q)
    kl = lambda a, b: float((a[a > 0] * np.log2(a[a > 0] / b[a > 0])).sum())
    return math.sqrt(max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m)))
def ent(p):
    p = p / p.sum(); p = p[p > 0]; return float(-(p * np.log2(p)).sum())
def drift(u, rg, tg):
    r, t = U[u][rg], U[u][tg]
    keys = sorted(set(r["c"]) | set(t["c"]))
    pr = np.array([r["c"].get(k, 0) for k in keys], float)
    pt = np.array([t["c"].get(k, 0) for k in keys], float)
    sr, st = set(r["c"]), set(t["c"])
    tr = max(r["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    tt = max(t["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    return dict(js_cells=js(pr, pt), dissim=1 - len(sr & st) / len(sr | st),
                anchor=float(tr != tt),
                novelty=sum(v for k, v in t["c"].items() if k not in sr) / sum(t["c"].values()),
                d_entropy=abs(ent(pt) - ent(pr)), js_slots=js(r["s"], t["s"]))
COMPS = ["js_cells", "dissim", "anchor", "novelty", "d_entropy", "js_slots"]
dec = {k: np.zeros(N) for k in COMPS}
for i, u in enumerate(users):
    a = drift(u, "RD", "TD")
    for k in COMPS: dec[k][i] = a[k]
cov = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)
comp = sum(np.argsort(np.argsort(dec[k])).astype(float) / (N - 1)
           for k in ["js_cells", "dissim", "anchor", "novelty"]) / 4

# ------------------------------------------------------------------ calibration
perm = rng.permutation(N); cal = perm[: N // 2]; ev = perm[N // 2:]
def nll(p, yv): return float(-np.mean(yv*np.log(np.clip(p,1e-12,1)) + (1-yv)*np.log(np.clip(1-p,1e-12,1))))
def brier(p, yv): return float(np.mean((p - yv) ** 2))
def ece(p, yv, nb=15):
    e = 0.0
    for k in range(nb):
        m = (p >= k/nb) & (p < (k+1)/nb if k < nb-1 else p <= 1.0)
        if m.sum(): e += m.sum()/len(p) * abs(p[m].mean() - yv[m].mean())
    return float(e)

soft = {t: np.exp(s1/t - lse[:, i]) for i, t in enumerate(TGRID)}
best = min(TGRID, key=lambda t: nll(soft[t][cal], y[cal]))
print(f"\nsoftmax temperature chosen on calibration users: T = {best}")
print("  T      NLL(cal)   Brier(cal)  mean p   accuracy")
for t in TGRID:
    print(f"  {t:<6} {nll(soft[t][cal],y[cal]):>8.4f} {brier(soft[t][cal],y[cal]):>11.4f} "
          f"{soft[t][cal].mean():>8.3f} {y[cal].mean():>10.3f}")
p_soft = soft[best]

def feats(k=K):
    return np.column_stack([s1, s1 - s2, s2 - topk[:, 9], s1 - topk[:, :k].mean(1), np.log(cov)])
X = feats()
mu, sd = X[cal].mean(0), X[cal].std(0)
Xz = (X - mu) / sd
def fit(Xf, yv, iters=80):
    Xf = np.column_stack([np.ones(len(yv)), Xf]); w = np.zeros(Xf.shape[1])
    for _ in range(iters):
        p = 1/(1+np.exp(-np.clip(Xf@w, -30, 30)))
        g = Xf.T@(p-yv)/len(yv)
        W = np.clip(p*(1-p), 1e-6, None)
        H = Xf.T@(Xf*W[:, None])/len(yv) + 1e-7*np.eye(Xf.shape[1])
        w = w - np.clip(np.linalg.solve(H, g), -5, 5)
    return w
w = fit(Xz[cal], y[cal])
p_platt = 1/(1+np.exp(-np.clip(np.column_stack([np.ones(N), Xz])@w, -30, 30)))
print(f"\nPlatt coefficients (s1, margin, s2-s10, s1-mean(top{K}), log coverage): {np.round(w,3)}")

print(f"\nEVALUATION users only (n = {len(ev):,}); accuracy = {y[ev].mean():.4f}")
print(f"  {'confidence':<12} {'NLL':>8} {'Brier':>8} {'ECE':>8} {'mean p':>8}")
VAR = {"margin (score)": (s1-s2)/np.maximum(s1,1e-12), f"softmax T={best}": p_soft, "platt": p_platt}
for n_, v in VAR.items():
    print(f"  {n_:<12} {nll(v[ev],y[ev]):>8.4f} {brier(v[ev],y[ev]):>8.4f} {ece(v[ev],y[ev]):>8.4f} {v[ev].mean():>8.3f}")

print(f"\nreliability of '{'platt'}' on evaluation users:")
print("  bin        n   mean p   accuracy   gap")
for k in range(10):
    m = (p_platt[ev] >= k/10) & (p_platt[ev] < (k+1)/10 if k < 9 else p_platt[ev] <= 1)
    if m.sum(): print(f"  {k/10:.1f}-{k/10+0.1:.1f} {m.sum():>6} {p_platt[ev][m].mean():>8.3f} "
                      f"{y[ev][m].mean():>10.3f} {y[ev][m].mean()-p_platt[ev][m].mean():>+7.3f}")

print("\nO1 PREVIEW -- calibration of the drift-blind calibrator, by drift quintile (eval users)")
q = np.quantile(comp[ev], [0, .2, .4, .6, .8, 1.0])
print(f"  {'quintile':<9} {'n':>6} {'mean drift':>11} {'accuracy':>9} {'mean p':>8} {'gap':>8} {'ECE':>7} {'Brier':>7}")
for k in range(5):
    lo, hi = q[k], q[k+1]
    m = (comp[ev] >= lo) & (comp[ev] <= hi if k == 4 else comp[ev] < hi)
    pe, ye = p_platt[ev][m], y[ev][m]
    print(f"  Q{k+1:<8} {m.sum():>6} {comp[ev][m].mean():>11.3f} {ye.mean():>9.3f} {pe.mean():>8.3f} "
          f"{ye.mean()-pe.mean():>+8.3f} {ece(pe,ye):>7.4f} {brier(pe,ye):>7.4f}")

def spear(a, b):
    """Tie-corrected Spearman (midranks). An argsort-of-argsort implementation is WRONG
    whenever either variable has ties -- the outcome y is binary, so it is all ties, and the
    error attenuated every r(., y) in the 2026-11-02 runs by 0.10 to 0.19."""
    from scipy import stats
    return float(stats.spearmanr(a, b).statistic)
print(f"\n  r(drift, p_platt) = {spear(comp,p_platt):+.3f}   r(drift, y) = {spear(comp,y):+.3f}   "
      f"r(drift, coverage) = {spear(comp,cov):+.3f}")
print(f"  truncation check: softmax over top-16 vs exact full pool at T={best}: "
      f"max abs diff {np.abs(np.exp(s1/best - (topk[:,:16].max(1)/best + np.log(np.exp(topk[:,:16]/best - topk[:,:16].max(1,keepdims=True)/best).sum(1)))) - p_soft).max():.4f}")

np.savez(f"results/cal_{LABEL}.npz", y=y, s1=s1, s2=s2, topk=topk, rank=rank, strue=strue,
         cov=cov, comp=comp, p_soft=p_soft, p_platt=p_platt, cal=cal, ev=ev, T=best, w=w,
         **{"dec_"+k: dec[k] for k in COMPS})
print(f"\nsaved results/cal_{LABEL}.npz")
