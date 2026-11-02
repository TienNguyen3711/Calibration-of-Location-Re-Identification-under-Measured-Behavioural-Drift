"""A second matcher family, and whether the mechanism reproduces under it.

WHY THIS EXISTS. `../../paper/OUTLINE.md` §6 names a single matcher family as the most likely
reason this paper is rejected, and C2b gives the objection a sharp form: if the mechanism
("measured drift carries information the matcher's own score lacks, but only once the
population has been disrupted") is real, it must not be an artefact of one profile
representation.

TWO ALTERNATIVE FAMILIES. Both hold the SCORING RULE fixed -- square-root transform, IDF
weighting, L2 normalisation, cosine -- and change only what the profile is made of, so any
difference attributes to the representation rather than to the similarity function.

  stbands   cell x six-hour band. The primary family profiles WHERE a person goes and throws
            away WHEN; this one keeps the joint. Requires a separate awk pass
            (`stbands_profiles.awk`). Genuinely more informative, and a different family.
  binary    cell presence, counts discarded. The primary family weights a cell by visit mass;
            this one records only that the cell was ever visited. Strictly less informative,
            and reuses the existing profiles -- it is the cheap end of the comparison.

WHAT IS TESTED, in both arms, for each family:
  - top-1 accuracy and the calibration of the Platt confidence (does the family work at all);
  - the reliability table by drift quintile (does the overconfidence collapse reappear);
  - the C2b increment: dAUC of drift beyond confidence and coverage, which the primary family
    puts at +0.026..+0.030 in the shock arm and approximately zero in the control.

PRE-REGISTERED, written before the run. PRIMARY ENDPOINT: the C2b increment per family per arm.

  REPRODUCES   both alternative families give a shock-arm increment whose bootstrap CI excludes
               zero, AND a control-arm increment materially smaller than it.
  FAILS        either family gives a shock-arm increment whose CI contains zero. The mechanism
               is representation-specific and the paper must say so.
  MIXED        anything else, reported per family rather than averaged.

Drift is NOT recomputed per family: it stays the primary cell-histogram Jensen-Shannon measure
on the odd-day half. That is deliberate. The question is whether a drift measure defined once
predicts failures of a DIFFERENT matcher, which is a stronger test than letting each family
have its own bespoke drift.

Usage: matcher_family2.py <windows.tsv> <stbands.tsv> <label> [pool] [seed]
"""
import sys, os, math, numpy as np
from scipy import sparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M

WPROF, SPROF, LABEL = sys.argv[1], sys.argv[2], sys.argv[3]
POOL = int(sys.argv[4]) if len(sys.argv) > 4 else 0
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 20261102
NB = 2000
REF, TGT = "W3", "W4"

print(f"=== second matcher family: {LABEL} ===")
U = M.parse_windows(WPROF, keep={REF+"M", REF+"D", TGT+"M", TGT+"D"})
S = M.parse_windows(SPROF, keep={REF+"M", TGT+"M"})
users = sorted([u for u, g in U.items()
                if len(g) == 4 and u in S and {REF+"M", TGT+"M"} <= set(S[u])], key=int)
if POOL and POOL < len(users):
    pick = np.random.default_rng(SEED).choice(len(users), POOL, replace=False)
    users = [users[i] for i in sorted(pick)]
N = len(users)
rng = np.random.default_rng(SEED)
cal = rng.permutation(N); ev = cal[N//2:]; cal = cal[:N//2]
print(f"pool {N:,}  cal {len(cal):,}  eval {len(ev):,}\n")

# drift and coverage, from the PRIMARY definition, on the odd-day half
d = M.drift(U, users, REF+"D", TGT+"D")
js = d["js_cells"]
cov = np.array([U[u][REF+"D"]["n"] + U[u][TGT+"D"]["n"] for u in users], float)

def binary_view(src):
    """Counts discarded: every visited cell gets weight 1."""
    return {u: {g: dict(n=v["n"], d=v["d"], n27=v.get("n27", 0),
                        c={k: 1 for k in v["c"]}, s=v["s"])
                for g, v in gs.items()} for u, gs in src.items()}

FAMILIES = {
    "primary (cells)":      (U, REF+"M", TGT+"M"),
    "stbands (cell x 6h)":  (S, REF+"M", TGT+"M"),
    "binary (presence)":    (binary_view(U), REF+"M", TGT+"M"),
}

def auc(p, yv):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = yv.sum(); n0 = len(yv) - n1
    return float((r[yv == 1].sum() - n1*(n1+1)/2)/(n1*n0))
def z(v): return (v - v.mean())/max(v.std(), 1e-12)
boot = rng.integers(0, len(ev), size=(NB, len(ev)))

OUT = {}
for fam, (src, rg, tg) in FAMILIES.items():
    nvoc = len({k for u in users for k in src[u][tg]["c"]})
    res = M.match(src, users, rg, tg)
    k = M.calibrate(res, cov, cal, ev)
    y, p = res["y"], k["p_platt"]
    base = np.column_stack([z(p), z(np.log(cov))])
    full = np.column_stack([base, z(js)])
    a = []
    for X in (base, full):
        w = M.logit_fit(X[cal], y[cal]); pp = M.logit_apply(X, w)
        a.append((auc(pp[ev], y[ev]), pp))
    bd = np.array([auc(a[1][1][ev][b], y[ev][b]) - auc(a[0][1][ev][b], y[ev][b]) for b in boot])
    lo, hi = np.quantile(bd, [0.025, 0.975])
    OUT[fam] = dict(y=y, p=p, top1=float(y.mean()), ece=float(M.ece(p[ev], y[ev])),
                    brier=float(M.brier(p[ev], y[ev])),
                    dauc=a[1][0]-a[0][0], lo=float(lo), hi=float(hi), nvoc=nvoc)
    print(f"{fam:<22} vocab {nvoc:>7,}  top-1 {y.mean():.4f}  "
          f"ECE {M.ece(p[ev],y[ev]):.4f}  Brier {M.brier(p[ev],y[ev]):.4f}")

print("\nreliability by DRIFT quintile (eval users), gap = accuracy - mean confidence")
q = np.quantile(js[ev], [0, .2, .4, .6, .8, 1.0])
masks = [(js[ev] >= q[i]) & (js[ev] <= q[i+1] if i == 4 else js[ev] < q[i+1]) for i in range(5)]
print(f"  {'family':<22} " + " ".join(f"{'Q'+str(i+1):>8}" for i in range(5)) + f" {'dcECE':>8}")
for fam, o in OUT.items():
    gaps = [o["y"][ev][m].mean() - o["p"][ev][m].mean() for m in masks]
    dce = sum(m.sum()/len(ev)*M.ece(o["p"][ev][m], o["y"][ev][m]) for m in masks)
    print(f"  {fam:<22} " + " ".join(f"{g:>+8.3f}" for g in gaps) + f" {dce:>8.4f}")

print("\nPRIMARY ENDPOINT: C2b increment, dAUC of drift beyond confidence and coverage")
print(f"  {'family':<22} {'dAUC':>8}  {'95% CI':>20}")
for fam, o in OUT.items():
    star = "" if o["lo"] <= 0 <= o["hi"] else "  CI excludes 0"
    print(f"  {fam:<22} {o['dauc']:>+8.4f}  [{o['lo']:>+7.4f},{o['hi']:>+7.4f}]{star}")

np.savez(f"results/family2_{LABEL}.npz",
         **{f"{k.split()[0]}_{f}": np.asarray(o[f])
            for k, o in OUT.items() for f in ("y", "p", "dauc", "lo", "hi")},
         js=js, cov=cov, cal=cal, ev=ev)
print(f"\nsaved results/family2_{LABEL}.npz")
