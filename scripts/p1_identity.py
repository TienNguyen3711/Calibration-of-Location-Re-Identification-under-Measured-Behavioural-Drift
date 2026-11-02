"""P1: bound the prevalence of pseudonym recycling, which A_PREMISES calls 'not verifiable'.

It is not fully verifiable, but the specific failure mode is bounded. If a uid were handed to
a different person mid-period, that user's self-similarity between consecutive windows would
show a sharp BREAK: high before, high after, one very low value across the switch. Continuous
behavioural drift instead produces a gradual decline with no isolated break.

Two quantities:
  (a) self-similarity across the four consecutive window pairs, against a null built from
      random user pairs in the same window pair;
  (b) the share of users showing a break pattern -- one pair below the null median while the
      others stay above a high threshold. That share is the UPPER BOUND on recycling, because
      ordinary drift also produces it sometimes.

Usage: p1_identity.py <windows.tsv> <label> [n_users] [seed]
"""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M
from scipy import sparse

PROF, LABEL = sys.argv[1], sys.argv[2]
NU = int(sys.argv[3]) if len(sys.argv) > 3 else 8000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261102

U = M.parse_windows(PROF)
need = {f"W{w}M" for w in range(5)}
comp = sorted([u for u, g in U.items() if need <= set(g)], key=int)
rng = np.random.default_rng(SEED)
users = [comp[i] for i in sorted(rng.choice(len(comp), min(NU, len(comp)), replace=False))]
N = len(users)
print(f"=== {LABEL} ===")
print(f"users with an M-half profile in all five windows: {len(comp):,}; using {N:,}\n")

vocab = {}
for u in users:
    for w in range(5):
        for k in U[u][f"W{w}M"]["c"]: vocab.setdefault(k, len(vocab))
V = len(vocab)
def build(w):
    rows, cols, vals = [], [], []
    for i, u in enumerate(users):
        c = U[u][f"W{w}M"]["c"]; tot = sum(c.values())
        for k, v in c.items():
            rows.append(i); cols.append(vocab[k]); vals.append(math.sqrt(v/tot))
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, V))
    n = np.sqrt(np.asarray(A.multiply(A).sum(1)).ravel()); n[n == 0] = 1
    return sparse.diags(1/n) @ A
W = [build(w) for w in range(5)]

PAIRS = [(0,1), (1,2), (2,3), (3,4)]
self_sim = np.zeros((N, len(PAIRS)))
null_med = np.zeros(len(PAIRS)); null_p95 = np.zeros(len(PAIRS))
for j, (a, b) in enumerate(PAIRS):
    self_sim[:, j] = np.asarray(W[a].multiply(W[b]).sum(1)).ravel()
    sh = rng.permutation(N)
    nullv = np.asarray(W[a].multiply(W[b][sh]).sum(1)).ravel()
    null_med[j] = np.median(nullv); null_p95[j] = np.quantile(nullv, .95)

print("Self-similarity of the SAME uid across consecutive windows, vs a shuffled-pair null")
print(f"  {'pair':<9} {'self median':>12} {'self p05':>9} {'null median':>12} {'null p95':>9} "
      f"{'self<null med':>14}")
for j, (a, b) in enumerate(PAIRS):
    s = self_sim[:, j]
    print(f"  W{a}->W{b}    {np.median(s):>12.4f} {np.quantile(s,.05):>9.4f} "
          f"{null_med[j]:>12.4f} {null_p95[j]:>9.4f} {100*(s<null_med[j]).mean():>13.2f}%")

print("\nShape of the decline: a recycled id gives an isolated break, drift a gradual slope.")
print(f"  the null is degenerate -- random pairs share almost no support, so the null MEDIAN is")
print(f"  {null_med.mean():.4f} and nothing can fall below it. The usable reference is the null p95")
print(f"  ({null_p95.mean():.4f} mean over the four pairs) and absolute levels.")
mn = self_sim.min(1)
print(f"\n  distribution of each user's WEAKEST consecutive-window self-similarity:")
print("    " + "  ".join(f"p{q:02d} {np.quantile(mn, q/100):.3f}" for q in (1, 5, 25, 50)))
print(f"\n  {'threshold':<34} {'users':>8} {'share':>8}")
for name, thr in (("null p95 (mean over pairs)", float(null_p95.mean())),
                  ("5x null p95", 5*float(null_p95.mean())),
                  ("absolute 0.05", 0.05), ("absolute 0.10", 0.10)):
    print(f"  weakest pair below {name:<16} {int((mn<thr).sum()):>8,} {100*(mn<thr).mean():>7.2f}%")
others_hi = np.array([np.sort(self_sim[i])[1:].min() for i in range(N)])
print(f"\n  BREAK pattern = weakest pair below 0.10 while the other three stay above a floor:")
for thr in (0.3, 0.4, 0.5):
    brk = (mn < 0.10) & (others_hi > thr)
    print(f"    other three > {thr}: {brk.sum():,} / {N:,} = {100*brk.mean():.3f}%")
print(f"  users whose EVERY pair is below 0.10: {int((self_sim<0.10).all(1).sum()):,} / {N:,} = "
      f"{100*(self_sim<0.10).all(1).mean():.3f}%")

pos = np.array([int(np.argmin(self_sim[i])) for i in range(N)])
print("\nWhere the weakest link falls (a mid-period re-pseudonymisation would cluster; the")
print("W3->W4 pair is the shock boundary, so mass there is expected and is NOT evidence of recycling):")
for j, (a, b) in enumerate(PAIRS):
    print(f"  weakest pair = W{a}->W{b}: {100*(pos==j).mean():>5.1f}%")
print("\nVERDICT: the bound above is an upper bound on pseudonym recycling. P1 stays an")
print("assumption; what changes is that its failure mode is now quantified, not unexamined.")
