"""Shared matcher / drift / calibration machinery.

Extracted so the sweeps (P3 gap, P6 pool, P7 coverage) run the SAME code as the primary
analysis. scripts/calibrate.py is deliberately left untouched as the record of the
2026-11-02 Arm A/B run; sweeps.py checks it reproduces calibrate.py at g = 0.
"""
import math, numpy as np
from scipy import sparse

def parse_windows(path, keep=None):
    """Read a windows_profiles.awk file. Returns {uid: {group: {...}}}."""
    U = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        g = f[2]
        if keep and g not in keep: continue
        cells = {}
        for tok in f[6][1:].split():
            k, v = tok.rsplit(":", 1); cells[k] = int(v)
        slots = np.zeros(48)
        for tok in f[7][1:].split():
            k, v = tok.rsplit(":", 1); slots[int(k)] = int(v)
        U.setdefault(f[1], {})[g] = dict(n=int(f[3]), d=int(f[4]), n27=int(f[5]),
                                         c=cells, s=slots)
    return U

# ----------------------------------------------------------------- matcher
def match(U, users, rg, tg, K=64, TGRID=(0.01,0.015,0.02,0.03,0.05,0.08,0.12,0.2,0.3), B=512):
    """Score every target profile in window `tg` against every reference profile in `rg`.
    Both are M-half groups. Returns top-K scores, true rank, and exact full-pool log-sum-exp
    on the temperature grid."""
    N = len(users)
    vocab = {}
    for u in users:
        for g in (rg, tg):
            for k in U[u][g]["c"]: vocab.setdefault(k, len(vocab))
    V = len(vocab)
    def build(grp):
        rows, cols, vals = [], [], []
        for i, u in enumerate(users):
            c = U[u][grp]["c"]; tot = sum(c.values())
            for k, v in c.items():
                rows.append(i); cols.append(vocab[k]); vals.append(math.sqrt(v/tot))
        return sparse.csr_matrix((vals, (rows, cols)), shape=(N, V))
    R, T = build(rg), build(tg)
    df = np.asarray((R > 0).sum(0)).ravel() + np.asarray((T > 0).sum(0)).ravel()
    Dg = sparse.diags(np.log((2*N + 1)/(df + 1)))
    R, T = R @ Dg, T @ Dg
    def l2(M):
        n = np.sqrt(np.asarray(M.multiply(M).sum(1)).ravel()); n[n == 0] = 1
        return sparse.diags(1/n) @ M
    R, T = l2(R), l2(T)
    Rt = R.T.tocsc()
    K = min(K, N - 1)
    topk = np.zeros((N, K)); rank = np.zeros(N, int); top = np.zeros(N, int)
    lse = np.zeros((N, len(TGRID)))
    for a in range(0, N, B):
        b = min(a + B, N)
        S = np.asarray((T[a:b] @ Rt).todense())
        idx = np.argpartition(-S, K, axis=1)[:, :K]
        vals = np.take_along_axis(S, idx, 1)
        order = np.argsort(-vals, axis=1)
        topk[a:b] = np.take_along_axis(vals, order, 1)
        top[a:b] = np.take_along_axis(idx, order, 1)[:, 0]
        st = S[np.arange(b - a), np.arange(a, b)]
        rank[a:b] = (S > st[:, None]).sum(1) + 1
        mx = S.max(1, keepdims=True)
        for ti, t in enumerate(TGRID):
            lse[a:b, ti] = (mx[:, 0]/t) + np.log(np.exp(S/t - mx/t).sum(1))
    y = (top == np.arange(N)).astype(float)
    return dict(y=y, topk=topk, rank=rank, lse=lse, TGRID=list(TGRID), K=K)

# ----------------------------------------------------------------- drift
def js(p, q):
    p, q = p/p.sum(), q/q.sum(); m = 0.5*(p + q)
    kl = lambda a, b: float((a[a > 0]*np.log2(a[a > 0]/b[a > 0])).sum())
    return math.sqrt(max(0.0, 0.5*kl(p, m) + 0.5*kl(q, m)))

def ent(p):
    p = p/p.sum(); p = p[p > 0]; return float(-(p*np.log2(p)).sum())

def drift_one(r, t):
    keys = sorted(set(r["c"]) | set(t["c"]))
    pr = np.array([r["c"].get(k, 0) for k in keys], float)
    pt = np.array([t["c"].get(k, 0) for k in keys], float)
    sr, st = set(r["c"]), set(t["c"])
    tr = max(r["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    tt = max(t["c"].items(), key=lambda kv: (kv[1], kv[0]))[0]
    return dict(js_cells=js(pr, pt), dissim=1 - len(sr & st)/len(sr | st),
                anchor=float(tr != tt),
                novelty=sum(v for k, v in t["c"].items() if k not in sr)/sum(t["c"].values()),
                d_entropy=abs(ent(pt) - ent(pr)), js_slots=js(r["s"], t["s"]))

COMPS = ["js_cells", "dissim", "anchor", "novelty", "d_entropy", "js_slots"]

def drift(U, users, rg, tg):
    out = {k: np.zeros(len(users)) for k in COMPS}
    for i, u in enumerate(users):
        d = drift_one(U[u][rg], U[u][tg])
        for k in COMPS: out[k][i] = d[k]
    return out

# ----------------------------------------------------------------- calibration
def nll(p, y):
    return float(-np.mean(y*np.log(np.clip(p, 1e-12, 1)) + (1-y)*np.log(np.clip(1-p, 1e-12, 1))))
def brier(p, y): return float(np.mean((p - y)**2))
def ece(p, y, nb=15):
    e = 0.0
    for k in range(nb):
        m = (p >= k/nb) & (p < (k+1)/nb if k < nb-1 else p <= 1.0)
        if m.sum(): e += m.sum()/len(p)*abs(p[m].mean() - y[m].mean())
    return float(e)
def spear(a, b):
    """Tie-corrected Spearman (midranks). An argsort-of-argsort implementation is WRONG
    whenever either variable has ties -- the outcome y is binary, so it is all ties, and the
    error attenuated every r(., y) in the 2026-11-02 runs by 0.10 to 0.19."""
    from scipy import stats
    return float(stats.spearmanr(a, b).statistic)

def platt_features(topk, cov, K=64):
    s1, s2 = topk[:, 0], topk[:, 1]
    k10 = min(9, topk.shape[1] - 1)
    return np.column_stack([s1, s1 - s2, s2 - topk[:, k10],
                            s1 - topk[:, :K].mean(1), np.log(cov)])

def logit_fit(X, y, iters=80):
    X = np.column_stack([np.ones(len(y)), X]); w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1/(1 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - y)/len(y)
        W = np.clip(p*(1 - p), 1e-6, None)
        H = X.T @ (X*W[:, None])/len(y) + 1e-7*np.eye(X.shape[1])
        w = w - np.clip(np.linalg.solve(H, g), -5, 5)
    return w
def logit_apply(X, w):
    return 1/(1 + np.exp(-np.clip(np.column_stack([np.ones(len(X)), X]) @ w, -30, 30)))

def calibrate(res, cov, cal, ev, K=64):
    """Fit softmax temperature and Platt on `cal`; return both probabilities for everyone."""
    s1 = res["topk"][:, 0]; y = res["y"]
    soft = {t: np.exp(s1/t - res["lse"][:, i]) for i, t in enumerate(res["TGRID"])}
    T = min(res["TGRID"], key=lambda t: nll(soft[t][cal], y[cal]))
    X = platt_features(res["topk"], cov, K)
    mu, sd = X[cal].mean(0), X[cal].std(0); sd[sd == 0] = 1
    Xz = (X - mu)/sd
    w = logit_fit(Xz[cal], y[cal])
    return dict(T=T, p_soft=soft[T], p_platt=logit_apply(Xz, w), w=w, mu=mu, sd=sd)
