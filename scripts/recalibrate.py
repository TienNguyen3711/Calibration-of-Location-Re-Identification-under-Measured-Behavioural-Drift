"""Drift-aware recalibration: the falsifiable prediction, run at last.

WHY THIS EXISTS. OUTLINE.md section 7.5 states a prediction that a reviewer can test and that
we had not run:

    drift-aware recalibration should improve calibration in the SHOCK arm and buy
    NOTHING in the control. If it helps in the control too, C2b is wrong.

C2b says the matcher's own score already absorbs drift in ordinary conditions, and that drift
only carries extra information once the score stops being a sufficient statistic for its own
correctness. If that is true, giving the calibrator access to drift must be useless in Arm B.

A second question, owed to `yuksekgonul2023beyond` and named as an open item in
results/ATYP_BASELINE.md section 7: their atypicality-aware recalibration (AAR) is a published
METHOD, and the atypicality baseline only tested whether rarity EXPLAINS the failure, never
whether their recalibrator FIXES it. Both are fitted here, in their own functional form, so the
comparison is to their method rather than to a strawman.

MODELS. All coefficients are fitted on CALIBRATION users only and evaluated on EVALUATION users
only, reusing the stored split so nothing is re-fitted that should not be.

  blind       the published Platt calibrator, read from the record. No drift, no atypicality.
  +drift      blind logit plus drift as a linear feature.
  AAR(atyp)   Yuksekgonul's form with atypicality: logit' = phi(a) * logit_blind + b,
              phi(a) = c2*a^2 + c1*a + c0, fitted by logistic regression on
              [logit, a*logit, a^2*logit]. An atypicality-dependent temperature.
  AAR(drift)  the same functional form with drift in place of atypicality. This is the
              drift-aware analogue of their recalibrator and is the candidate method.
  AAR(both)   phi(a, d) with both modulators, to see whether they are complementary.

PRE-REGISTERED, written before the run:

  PRIMARY ENDPOINT: the reduction in |gap| in the WORST drift quintile on evaluation users,
  AAR(drift) against blind, reported per arm with a 2,000-resample bootstrap CI.

  Reading (a), C2b SUPPORTED: Arm A's reduction CI excludes zero AND Arm B's reduction is less
  than half of Arm A's point estimate.
  Reading (b), C2b CONTRADICTED: Arm B's reduction is comparable to Arm A's (at least half the
  point estimate with a CI excluding zero). The paper's mechanism claim would have to be
  withdrawn, and this file is the record that we looked.
  Reading (c), method fails: Arm A's reduction CI contains zero. Drift predicts the failure but
  cannot be used to repair it, which is a weaker but publishable result and must be stated as
  such rather than quietly dropped.

  SECONDARY: aggregate ECE, drift-conditional ECE, Brier, NLL, AUC; and whether AAR(atyp)
  recovers any of what AAR(drift) recovers.

Nothing here can change the matcher or the drift measure; it only re-maps confidence.

Usage: recalibrate.py <results/cal_LABEL.npz> <results/atyp_LABEL.npz> <label>
"""
import sys, numpy as np
sys.path.insert(0, "scripts")
import matchlib as M

CAL, ATY, LABEL = sys.argv[1], sys.argv[2], sys.argv[3]
NB = 2000
SEED = 20261102

Z = np.load(CAL); A = np.load(ATY)
y, p_blind, cov = Z["y"], Z["p_platt"], Z["cov"]
js, comp, cal, ev = Z["dec_js_cells"], Z["comp"], Z["cal"], Z["ev"]
a_pop, a_knn, a_maha = A["a_pop"], A["a_knn"], A["a_maha"]
N = len(y)
rng = np.random.default_rng(SEED)
boot = rng.integers(0, len(ev), size=(NB, len(ev)))

print(f"=== drift-aware recalibration: {LABEL} ===")
print(f"pool {N:,}  cal {len(cal):,}  eval {len(ev):,}  accuracy {y.mean():.4f}")
print(f"blind calibrator on eval: ECE {M.ece(p_blind[ev],y[ev]):.4f}  "
      f"Brier {M.brier(p_blind[ev],y[ev]):.4f}  NLL {M.nll(p_blind[ev],y[ev]):.4f}\n")

def z(v): return (v - v[cal].mean()) / max(v[cal].std(), 1e-12)
logit = np.log(np.clip(p_blind, 1e-6, 1-1e-6) / (1 - np.clip(p_blind, 1e-6, 1-1e-6)))

def fit(X):
    """Fit on calibration users, return probabilities for everyone."""
    mu, sd = X[cal].mean(0), X[cal].std(0); sd[sd == 0] = 1
    Xz = (X - mu) / sd
    w = M.logit_fit(Xz[cal], y[cal])
    return M.logit_apply(Xz, w)

zd, za = z(js), z(a_pop)
zk, zm = z(a_knn), z(a_maha)
MODELS = {
    "blind (published)": p_blind,
    "+drift (linear)":   fit(np.column_stack([logit, zd])),
    "AAR(atyp)":         fit(np.column_stack([logit, za*logit, za**2*logit,
                                              zk*logit, zm*logit])),
    "AAR(drift)":        fit(np.column_stack([logit, zd*logit, zd**2*logit])),
    "AAR(both)":         fit(np.column_stack([logit, zd*logit, zd**2*logit,
                                              za*logit, za**2*logit])),
}

# ------------------------------------------------------------------ aggregate
print("EVALUATION users, aggregate")
print(f"  {'model':<20} {'ECE':>8} {'Brier':>8} {'NLL':>8} {'AUC':>8}")
def auc(p, yv):
    r = np.argsort(np.argsort(p)).astype(float) + 1
    n1 = yv.sum(); n0 = len(yv) - n1
    return float((r[yv == 1].sum() - n1*(n1+1)/2)/(n1*n0))
for k, p in MODELS.items():
    print(f"  {k:<20} {M.ece(p[ev],y[ev]):>8.4f} {M.brier(p[ev],y[ev]):>8.4f} "
          f"{M.nll(p[ev],y[ev]):>8.4f} {auc(p[ev],y[ev]):>8.4f}")

# ------------------------------------------------------------------ by drift quintile
q = np.quantile(js[ev], [0, .2, .4, .6, .8, 1.0])
masks = []
for k in range(5):
    masks.append((js[ev] >= q[k]) & (js[ev] <= q[k+1] if k == 4 else js[ev] < q[k+1]))

print("\nEVALUATION users, calibration gap (accuracy - mean confidence) by DRIFT quintile")
print(f"  {'model':<20} " + " ".join(f"{'Q'+str(i+1):>8}" for i in range(5)) + f" {'dcECE':>8}")
for k, p in MODELS.items():
    gaps = [y[ev][m].mean() - p[ev][m].mean() for m in masks]
    dce = sum(m.sum()/len(ev)*M.ece(p[ev][m], y[ev][m]) for m in masks)
    print(f"  {k:<20} " + " ".join(f"{g:>+8.3f}" for g in gaps) + f" {dce:>8.4f}")

# ------------------------------------------------------------------ PRIMARY
print("\nPRIMARY ENDPOINT: reduction in |gap| in the WORST drift quintile (Q5), vs blind")
worst = masks[4]
def absgap(p, idx):
    return abs(y[ev][idx].mean() - p[ev][idx].mean())
wi = np.where(worst)[0]
bw = rng.integers(0, len(wi), size=(NB, len(wi)))
base_g = absgap(p_blind, worst)
print(f"  blind |gap| in Q5 = {base_g:.4f}  (n = {worst.sum():,})")
print(f"  {'model':<20} {'|gap|':>8} {'reduction':>10} {'95% CI':>20}")
RES = {}
for k, p in MODELS.items():
    if k == "blind (published)": continue
    g = absgap(p, worst)
    red = base_g - g
    bs = np.array([abs(y[ev][wi[b]].mean() - p_blind[ev][wi[b]].mean())
                   - abs(y[ev][wi[b]].mean() - p[ev][wi[b]].mean()) for b in bw])
    lo, hi = np.quantile(bs, [0.025, 0.975])
    RES[k] = (red, lo, hi)
    star = "" if lo <= 0 <= hi else "  CI excludes 0"
    print(f"  {k:<20} {g:>8.4f} {red:>+10.4f}  [{lo:>+7.4f},{hi:>+7.4f}]{star}")

np.savez(f"results/recal_{LABEL}.npz",
         **{("p_" + k.replace(" ", "_").replace("(", "").replace(")", "")): v
            for k, v in MODELS.items()},
         **{("red_" + k.replace(" ", "_").replace("(", "").replace(")", "")): np.array(v)
            for k, v in RES.items()})
print(f"\nsaved results/recal_{LABEL}.npz")
