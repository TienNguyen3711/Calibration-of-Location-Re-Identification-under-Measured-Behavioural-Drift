# The atypicality baseline

Run 2026-11-02. Discharges the obligation created by
[`R1_ATYPICALITY_RULING.md`](../../../Material/_ai_research/intake/R1_ATYPICALITY_RULING.md) §6,
which was itself created by reading Yuksekgonul, Zhang, Zou & Guestrin, *Beyond Confidence:
Reliable Models Should Also Consider Atypicality*, NeurIPS 2023 (`10.52202/075280-1669`) in full.

**Pre-registered verdict: reading (a), in its strongest form. Drift is not explained by rarity.**
The primary endpoint is **ΔAUC(drift | base + atypicality) = +0.0560 [+0.0509, +0.0612]** in the
shock arm — *larger* than the published C2b increment over confidence alone (+0.0307) — against
**+0.0043 [+0.0002, +0.0087]** in the matched control. Atypicality on its own adds **nothing** in
either arm. One honest complication, which is a finding rather than a problem, is in §5.

## 1. Why this run exists

Their Table 1 already reports ECE **by per-sample atypicality quantile**, where

    a_X(x) = − max_y log P(X = x | Y = y),   "with respect to the training distribution".

So the *shape* of our O1 §4b analysis is theirs. The question a reviewer who knows that paper
asks is whether a plain rarity score of the **target window alone** — one window, no reference
window, no parity split — explains the calibration failure just as well. If it does, `δ` buys
nothing and the paper has no contribution.

Pre-registered before the run, in the ruling §6: (a) both predict and `δ` retains an independent
coefficient → contribution stands; (b) atypicality absorbs `δ` → the paper shrinks to
"Yuksekgonul transfers to privacy attacks"; (c) `δ` absorbs atypicality → strongest outcome.
Decision rule and primary endpoint were fixed in the script header before execution
(`../scripts/atyp_baseline.py`).

## 2. Method

Nothing was refitted. The matcher, the Platt calibrator, the cal/eval split and `δ` are read from
`cal_armA_dataset2.npz` and `cal_armB_dataset1.npz`, the record of the 2026-11-02 Arm A/B run.
User order is reconstructed from the profile file by `calibrate.py`'s own rule, including Arm B's
subsample of 24,517 from 99,549 at seed 20261102, and then **checked against the stored coverage
and `js_cells` arrays: both matched at `max|Δ| = 0.000e+00` in both arms.** A wrong draw would
have shown up here.

Three estimators, spanning their menu, all computed on the **target window, matcher half (TM)** —
i.e. on the records the matcher *does* consume, which is the whole point of the baseline — with
the "training distribution" being the population of reference-window matcher halves (RM):

| | What it is | Needs the matcher? |
|---|---|---|
| `a_pop` | mean negative log population probability of the user's visited cells; their Def. 2.1 collapsed to a single class, `−log P(x)` | **No.** No scores, no embedding, no model |
| `a_knn` | their KNN estimator: `1 − mean(top-5 cosine)`, read off the stored top-64 scores | Yes — it *is* the matcher's neighbourhood |
| `a_maha` | their vision estimator: one Gaussian, tied covariance, on the predictor's representation reduced by truncated SVD; Mahalanobis² from the reference-population mean | Yes — the representation |

`base = [z(p_platt), z(log coverage)]`, identical to the C2b nested fit in `sweeps.py`, so the
numbers are directly comparable. Coefficients are fitted on the calibration half; AUC and the
2,000-resample bootstrap are on the evaluation half only.

Reproduction check: this run gives `drift | base` = **+0.0307** where the published per-gap sweep
gives **+0.0292** at g = 0, with NLL 0.4209 → 0.3731 against the published 0.4203 → 0.3733. The
small AUC difference is the cal/eval split (`sweeps.py` draws its own); the mechanism reproduces.

## 3. The decisive table: ΔAUC in both arms

| Increment | **Arm A (shock)** | **Arm B (control)** |
|---|---|---|
| atypicality (all three) \| base | **−0.0008** [−0.0016, −0.0000] | **−0.0030** [−0.0059, −0.0004] |
| `a_pop` \| base (the model-free one) | −0.0011 [−0.0016, −0.0005] | −0.0031 [−0.0060, −0.0006] |
| `a_knn` \| base | −0.0006 [−0.0011, −0.0001] | +0.0002 [−0.0001, +0.0004] |
| `a_maha` \| base | +0.0002 [+0.0001, +0.0003] | −0.0013 [−0.0023, −0.0004] |
| drift \| base *(published C2b)* | +0.0307 [+0.0264, +0.0348] | −0.0045 [−0.0083, −0.0008] |
| **drift \| base + atypicality** ← primary | **+0.0560 [+0.0509, +0.0612]** | **+0.0043 [+0.0002, +0.0087]** |
| atypicality \| base + drift | +0.0245 [+0.0212, +0.0278] | +0.0058 [+0.0026, +0.0091] |

Every rarity measure, including the two that need the matcher, adds **nothing or slightly less
than nothing** over confidence and coverage. The small negative values are the out-of-sample cost
of a useless parameter, not a signal.

Drift's increment **survives conditioning on all three rarity measures and roughly doubles**
(+0.0307 → +0.0560). And the **arm contrast is preserved at a factor of 13**: +0.0560 under the
shock against +0.0043 in the control.

**This refutes reading (b).** Rarity does not absorb drift.

Dimension sensitivity, in fairness to their estimator — `a_maha` at d = 64 captures only 7.5% of
the reference energy, so it was given 4× and 8× more (`atyp_armA_svd_sensitivity.txt`):

| d | SVD energy | `a_maha` \| base | primary endpoint |
|---|---|---|---|
| 64 | 0.075 | +0.0002 [+0.0001, +0.0003] | +0.0560 |
| 256 | 0.165 | +0.0008 [+0.0003, +0.0012] | +0.0559 |
| 512 | 0.242 | +0.0005 [+0.0001, +0.0009] | +0.0561 |

The estimator was not crippled by its dimension; eight times the dimensions changes nothing.

## 4. The reliability tables: the collapse lives on the drift axis and nowhere else

Their Table 1 shape, eval users of the shock arm. `gap` = accuracy − mean confidence, so negative
means **overconfident**:

| Stratifier | Q1 gap | Q5 gap | Q1 ECE | Q5 ECE |
|---|---|---|---|---|
| `a_pop` (model-free rarity) | −0.002 | +0.006 | 0.0210 | 0.0208 |
| `a_knn` | +0.014 | +0.019 | 0.0233 | 0.0403 |
| `a_maha` | −0.028 | +0.010 | 0.0332 | 0.0341 |
| **`js_cells` (drift, ours)** | **+0.044** | **−0.202** | **0.0444** | **0.2018** |
| `comp` (composite drift) | +0.054 | −0.175 | 0.0544 | 0.1752 |

The full drift ladder is monotone and ends in collapse: gap +0.044, +0.069, +0.073, +0.041,
**−0.202**; ECE 0.044, 0.069, 0.073, 0.067, **0.2018**. No rarity axis has anything resembling
it — the largest gap any of them produces is **+0.019**, and it is *under*-confidence.

In the control arm the drift ladder flattens exactly as C2b predicts: gap −0.001, +0.009, +0.012,
+0.013, −0.010, with worst-quintile ECE 0.0225 against 0.2018 under the shock.

So: **the overconfidence collapse appears on a drift axis, in the shock arm, and nowhere else.**

## 5. The honest complication, which goes in the paper

`atypicality | base + drift` is positive in both arms (+0.0245, +0.0058), and the full model beats
base + drift on eval users in Arm A (AUC 0.9130 against 0.8885, NLL 0.3320 against 0.3731). That
is a genuine out-of-sample gain, with coefficients fitted on the calibration half, not overfitting.

It is driven by `a_knn`, and the mechanism is plain in the residual regression (T3, Arm A):

| | alone | jointly with `js_cells` |
|---|---|---|
| `a_knn` | +0.0015 [−0.0050, +0.0084] | **+0.1160** [+0.1060, +0.1259] |
| `js_cells` | −0.0746 [−0.0800, −0.0688] | **−0.1550** [−0.1636, −0.1463] |

`a_knn` correlates **+0.709** with `js_cells` (+0.741 in Arm B) and is null on its own. Allowed
into the same model with opposite signs, the collinear pair separates *how good the match is in
absolute terms* from *how much the person changed*, and the contrast predicts better than either.

**What this is not.** It is not evidence that rarity explains the calibration failure: `a_knn` is
read off the matcher's own top-5 cosines, so it is not an exogenous covariate, and the one
estimator that is exogenous — `a_pop` — contributes −0.0011 in the shock arm and −0.0031 in the
control. The paper must report the suppression effect, state that it comes from a score-space
quantity, and not dress it up as an independent rarity signal.

## 6. A clean secondary fact worth a sentence in §2

Rarity and drift are **different axes with opposite signs**. `r(a_pop, y) = +0.082` in the shock
arm and +0.088 in the control: people who go to unusual places are *easier* to re-identify, as
distinctiveness should make them. `r(js_cells, y) = −0.537` and −0.315: people who changed are
*harder*. And the two are nearly orthogonal, `r(a_pop, js_cells) = +0.003` (A) and +0.036 (B).

A measure of how unusual someone is tells you almost nothing about how much they changed. That is
the one-sentence statement of why the covariate had to be a change measure.

## 7. What this does not establish

- **Their analysis shape remains theirs.** Nothing here makes stratified reliability analysis
  ours; `yuksekgonul2023beyond` is cited as the method antecedent and §2 must say so.
- **No atypicality-aware recalibration was fitted.** We tested whether rarity *explains* the
  failure, not whether their AAR *fixes* it. Their recalibrator conditions on `a(X)`; a
  drift-aware recalibrator is the falsifiable prediction already stated in `OUTLINE.md` §7.5, and
  neither has been run.
- **`a_pop` is one instantiation.** A different population density — a mixture, a POI-weighted
  prior, a per-slot model — could behave differently. Three estimators spanning their menu all
  giving ≈0 is strong, not conclusive.
- **Arm C is untouched.** At n = 56, calibration is not assessable there, so this baseline exists
  for YJMob100K only, which is the same limitation C2/C2b/C3 already carry (gap G6).
- **The run is on the target window's matcher half.** A rarity score computed on the drift half
  would be a different, and more favourable-to-us, comparison; it was deliberately not used,
  because the point of the baseline is the cheap alternative that needs no parity split.

## 8. Reproduce

```
gzcat data/raw/yjmob100k/yjmob100k-dataset2.csv.gz | awk -f scripts/p2_profiles.awk > prof_d2.tsv
gzcat data/raw/yjmob100k/yjmob100k-dataset1.csv.gz | awk -f scripts/p2_profiles.awk > prof_d1.tsv
python3 scripts/atyp_baseline.py prof_d2.tsv results/cal_armA_dataset2.npz armA_dataset2
python3 scripts/atyp_baseline.py prof_d1.tsv results/cal_armB_dataset1.npz armB_dataset1 24517
for d in 256 512; do python3 scripts/atyp_baseline.py prof_d2.tsv \
    results/cal_armA_dataset2.npz armA_d$d 0 $d; done
```

Outputs: `atyp_armA_output.txt`, `atyp_armB_output.txt`, `atyp_armA_svd_sensitivity.txt`,
`atyp_arm{A,B}_*.npz`. The profile `.tsv` files are regenerated rather than kept; the two awk
passes read 29.4M and 111.5M rows and take a few minutes each.

Note on the console: `matchlib.py` emits spurious `RuntimeWarning: divide by zero / overflow /
invalid value encountered in matmul` from the Accelerate BLAS on this machine. Every feature was
checked finite (`p_platt`, `log cov`, `js_cells`, `a_pop`, `a_knn`, `a_maha`) and the C2b
reproduction in §2 lands where it should, so the warnings are an artefact, not a defect.
