"""P2: is measured drift mechanically coupled to the matcher score?

Design under test (A_PREMISES P2):
  interleaved-day split. Within every window, even day indices (half M) feed the
  MATCHER, odd day indices (half D) feed the DRIFT measure. No record is used by
  both. The pre-registered diagnostic computes drift BOTH ways -- on the held-out
  half (decoupled, the design) and on the matcher's own half (coupled, the
  circular counterfactual) -- and compares their correlation with the score.
"""
import sys, math, numpy as np
from scipy import sparse

PROF = sys.argv[1]
OUT  = sys.argv[2] if len(sys.argv) > 2 else "p2_vars.npz"

def parse(path):
    U = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        uid, grp, nrec, nday = f[1], f[2], int(f[3]), int(f[4])
        cells = {}
        for tok in f[5][1:].split():
            k, v = tok.rsplit(":", 1); cells[k] = int(v)
        slots = np.zeros(48)
        for tok in f[6][1:].split():
            k, v = tok.rsplit(":", 1); slots[int(k)] = int(v)
        U.setdefault(uid, {})[grp] = dict(n=nrec, d=nday, c=cells, s=slots)
    return U

U = parse(PROF)
users = sorted([u for u, g in U.items() if {"RM","RD","TM","TD"} <= set(g)], key=int)
N = len(users)
print(f"users with all four (window x half) groups: {N:,} of {len(U):,}\n")

# ---------------------------------------------------------------- matcher (M half only)
vocab = {}
for u in users:
    for g in ("RM", "TM"):
        for k in U[u][g]["c"]:
            vocab.setdefault(k, len(vocab))
V = len(vocab)

def build(grp):
    rows, cols, vals = [], [], []
    for i, u in enumerate(users):
        c = U[u][grp]["c"]; tot = sum(c.values())
        for k, v in c.items():
            rows.append(i); cols.append(vocab[k]); vals.append(math.sqrt(v / tot))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(N, V))

R = build("RM"); T = build("TM")
df = np.asarray((R > 0).sum(axis=0)).ravel() + np.asarray((T > 0).sum(axis=0)).ravel()
idf = np.log((2 * N + 1) / (df + 1))
D = sparse.diags(idf)
R = R @ D; T = T @ D
def l2(M):
    n = np.sqrt(np.asarray(M.multiply(M).sum(axis=1)).ravel()); n[n == 0] = 1
    return sparse.diags(1 / n) @ M
R = l2(R); T = l2(T)

Rt = R.T.tocsc()
s1 = np.zeros(N); s2 = np.zeros(N); strue = np.zeros(N); rank = np.zeros(N, int); top = np.zeros(N, int)
B = 512
for a in range(0, N, B):
    b = min(a + B, N)
    S = np.asarray((T[a:b] @ Rt).todense())
    for j in range(b - a):
        row = S[j]; i = a + j
        o = np.argpartition(-row, 2)[:2]; o = o[np.argsort(-row[o])]
        s1[i], s2[i], top[i] = row[o[0]], row[o[1]], o[0]
        strue[i] = row[i]
        rank[i] = int((row > row[i]).sum()) + 1
y = (top == np.arange(N)).astype(float)
c_margin = (s1 - s2) / np.maximum(s1, 1e-12)
print(f"matcher (M half, cosine on IDF-weighted Hellinger cell profiles, pool = {N:,})")
print(f"  top-1 accuracy {y.mean():.4f}   top-10 {(rank<=10).mean():.4f}   median rank of true user {np.median(rank):.0f}")
print(f"  score s1 mean {s1.mean():.3f}   margin confidence mean {c_margin.mean():.3f}\n")

# ---------------------------------------------------------------- drift components
def js(p, q):
    p = p / p.sum(); q = q / q.sum(); m = 0.5 * (p + q)
    def kl(a, b):
        mask = a > 0
        return float((a[mask] * np.log2(a[mask] / b[mask])).sum())
    return math.sqrt(max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m)))

def ent(p):
    p = p / p.sum(); p = p[p > 0]
    return float(-(p * np.log2(p)).sum())

def drift(u, rg, tg):
    r, t = U[u][rg], U[u][tg]
    keys = sorted(set(r["c"]) | set(t["c"]))
    pr = np.array([r["c"].get(k, 0) for k in keys], float)
    pt = np.array([t["c"].get(k, 0) for k in keys], float)
    sr, st = set(r["c"]), set(t["c"])
    inter = len(sr & st); uni = len(sr | st)
    topr = max(r["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    topt = max(t["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    nov = sum(v for k, v in t["c"].items() if k not in sr) / sum(t["c"].values())
    return dict(js_cells=js(pr, pt), dissim=1 - inter / uni, anchor=float(topr != topt),
                novelty=nov, d_entropy=abs(ent(pt) - ent(pr)), js_slots=js(r["s"], t["s"]))

COMPS = ["js_cells", "dissim", "anchor", "novelty", "d_entropy", "js_slots"]
dec = {k: np.zeros(N) for k in COMPS}
cou = {k: np.zeros(N) for k in COMPS}
for i, u in enumerate(users):
    a = drift(u, "RD", "TD")
    b = drift(u, "RM", "TM")
    for k in COMPS:
        dec[k][i] = a[k]; cou[k][i] = b[k]

cov = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)
cov_m = np.array([U[u]["RM"]["n"] + U[u]["TM"]["n"] for u in users], float)

def spear(a, b):
    """Tie-corrected Spearman (midranks). An argsort-of-argsort implementation is WRONG
    whenever either variable has ties -- the outcome y is binary, so it is all ties, and the
    error attenuated every r(., y) in the 2026-11-02 runs by 0.10 to 0.19."""
    from scipy import stats
    return float(stats.spearmanr(a, b).statistic)

print("Spearman correlations (n = %d). c = margin confidence, y = top-1 correct.\n" % N)
hdr = f"{'component':<11} {'r(d_dec,c)':>11} {'r(d_cou,c)':>11} {'r(d_dec,y)':>11} {'r(d_cou,y)':>11} {'r(dec,cou)':>11} {'r(d_dec,cov)':>13}"
print(hdr); print("-" * len(hdr))
for k in COMPS:
    print(f"{k:<11} {spear(dec[k],c_margin):>11.3f} {spear(cou[k],c_margin):>11.3f} "
          f"{spear(dec[k],y):>11.3f} {spear(cou[k],y):>11.3f} {spear(dec[k],cou[k]):>11.3f} "
          f"{spear(dec[k],cov):>13.3f}")

def composite(d):
    z = np.zeros(N)
    for k in ["js_cells", "dissim", "anchor", "novelty", "js_slots"]:
        r = np.argsort(np.argsort(d[k])).astype(float) / (N - 1)
        z += r
    return z / 5
cd, cc = composite(dec), composite(cou)
print(f"\n{'composite':<11} {spear(cd,c_margin):>11.3f} {spear(cc,c_margin):>11.3f} "
      f"{spear(cd,y):>11.3f} {spear(cc,y):>11.3f} {spear(cd,cc):>11.3f} {spear(cd,cov):>13.3f}")
print(f"\ncoverage check: r(coverage_D, c) = {spear(cov,c_margin):.3f}, "
      f"r(coverage_M, c) = {spear(cov_m,c_margin):.3f}, r(cov_D, cov_M) = {spear(cov,cov_m):.3f}")

np.savez(OUT, y=y, c=c_margin, s1=s1, s2=s2, rank=rank, cov=cov, cov_m=cov_m,
         cd=cd, cc=cc, **{"dec_"+k: dec[k] for k in COMPS}, **{"cou_"+k: cou[k] for k in COMPS})
print(f"\nsaved {OUT}")

# ---------------------------------------------------------------- P2 follow-ups
print("\n" + "=" * 86)
print("Placebo: each user keeps their own score but is given ANOTHER user's drift")
print("=" * 86)
rng = np.random.default_rng(20261102)
for name, v in (("js_cells (decoupled)", dec["js_cells"]), ("composite (decoupled)", cd)):
    rc = [spear(v[rng.permutation(N)], c_margin) for _ in range(20)]
    ry = [spear(v[rng.permutation(N)], y) for _ in range(20)]
    print(f"  {name:<24} r(d,c) real {spear(v,c_margin):+.3f}  shuffled {np.mean(rc):+.4f} "
          f"(sd {np.std(rc):.4f}, max|.| {max(abs(x) for x in rc):.4f})")
    print(f"  {'':<24} r(d,y) real {spear(v,y):+.3f}  shuffled {np.mean(ry):+.4f} (sd {np.std(ry):.4f})")

print("\n" + "=" * 86)
print("Preliminary nested fit for y (this is O2; shown here only as a smell test)")
print("=" * 86)
def fit(X, yv, iters=60):
    X = np.column_stack([np.ones(len(yv)), X]); w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - yv) / len(yv)
        W = np.clip(p * (1 - p), 1e-6, None)
        H = X.T @ (X * W[:, None]) / len(yv) + 1e-7 * np.eye(X.shape[1])
        w = w - np.clip(np.linalg.solve(H, g), -5, 5)
    p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30)))
    return w, float(np.mean(yv*np.log(p+1e-12) + (1-yv)*np.log(1-p+1e-12))), p
def auc(p, yv):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = yv.sum(); n0 = len(yv) - n1
    return float((r[yv == 1].sum() - n1*(n1+1)/2) / (n1*n0))
def zs(v): return (v - v.mean()) / v.std()
lc = np.log(cov)
print(f"  {'model':<32} {'log-lik':>9} {'AUC':>7} {'Brier':>7}   coefficients")
for name, X in (("c alone", np.column_stack([zs(c_margin)])),
                ("c + coverage", np.column_stack([zs(c_margin), zs(lc)])),
                ("c + coverage + drift DECOUPLED", np.column_stack([zs(c_margin), zs(lc), zs(cd)])),
                ("c + coverage + drift COUPLED", np.column_stack([zs(c_margin), zs(lc), zs(cc)]))):
    w, ll, p = fit(X, y)
    print(f"  {name:<32} {ll:>+9.4f} {auc(p,y):>7.4f} {np.mean((p-y)**2):>7.4f}   {np.round(w,3)}")

print("\n  reliability of the drift measure (two disjoint day sets, same person):")
print(f"    test-retest r(js_cells odd days, js_cells even days) = {spear(dec['js_cells'],cou['js_cells']):.3f}")
print("\n  NOTE: c is a margin SCORE, not a probability (mean %.3f vs accuracy %.3f)." % (c_margin.mean(), y.mean()))
print("  Its raw Brier is therefore not a miscalibration measure. O1 needs a probability output.")
print("\n  reliability-diagram deciles of the raw margin score:")
b = np.clip((c_margin * 10).astype(int), 0, 9)
for k in range(10):
    m = b == k
    if m.sum(): print(f"    {k/10:.1f}-{k/10+0.1:.1f} | n={m.sum():>5} | mean c {c_margin[m].mean():.3f} | accuracy {y[m].mean():.3f}")
