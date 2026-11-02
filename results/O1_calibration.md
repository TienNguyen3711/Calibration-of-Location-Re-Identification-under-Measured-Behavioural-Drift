# O1: the matcher's probability output, and whether calibration survives drift

Run 2026-11-02. Covers two steps: giving `M` a probability output (which O1 was blocked on),
and rerunning the whole check on Arm B, the no-shock control.

Scripts: [../scripts/calibrate.py](../scripts/calibrate.py),
[../scripts/compare_arms.py](../scripts/compare_arms.py). Console output:
[O1_calibration_output.txt](O1_calibration_output.txt). Variables: `cal_armA_dataset2.npz`,
`cal_armB_dataset1.npz`.


> **Correction, 2026-11-02.** The first pass computed Spearman with an argsort-of-argsort
> implementation that ignores tied ranks. The outcome `y` is binary and the `anchor` component is
> binary, so every correlation involving either was attenuated — by 0.10 to 0.19 for `r(·, y)`. All
> such numbers below are the corrected, midrank values (`scipy.stats.spearmanr`). Correlations
> between two continuous variables were unaffected, so every ECE, calibration gap, ΔAUC and
> bootstrap interval in these documents is unchanged. `scripts/matchlib.py`, `p2_decoupling.py`,
> `calibrate.py` and `compare_arms.py` were all fixed and the affected analyses re-run.

Reproduce (from `Codebase/code/`, after building the profile files with `scripts/p2_profiles.awk`):
```
python3 scripts/calibrate.py prof_d2.tsv armA_dataset2
python3 scripts/calibrate.py prof_d1.tsv armB_dataset1 24517
python3 scripts/compare_arms.py
```

**Headline.** The drift–score association is **not** specific to the disruption: it is nearly as
strong in the control. What the shock breaks is **calibration**, and an aggregate calibration metric
hides that failure almost entirely. The interpretation in A_PREMISES §0 is narrowed accordingly
(§5 below).

## 1. Probability output

Two confidences are produced; neither sees any drift feature, by design — O1 asks whether a
drift-blind calibrator stays calibrated across drift strata.

- **softmax(T)**: `P = exp(s₁/T) / Σ_v exp(s_v/T)`, the exact closed-world posterior over the pool
  under a Gibbs model of the score. One parameter, chosen on a grid by NLL on calibration users.
- **platt**: logistic regression on score order statistics `s₁`, `s₁−s₂`, `s₂−s₁₀`,
  `s₁−mean(top 64)` and `log(coverage)`, fitted on calibration users.

Users are split 50/50 into calibration and evaluation sets. The **candidate pool is the full pool
in both cases** — the split decides whose metrics are reported, not who can be matched.

**The full-pool denominator is not optional.** Computing the softmax over the top 16 candidates
instead of all 24,517 differs from the exact value by up to **0.76** in probability at the chosen
`T = 0.03`. Any top-K shortcut here is wrong; `scripts/calibrate.py` accumulates the exact
log-sum-exp over the whole pool on a grid of `T`.

Arm A, evaluation users (n = 12,259), accuracy 0.7149:

| confidence | NLL | Brier | ECE | mean p |
|---|---|---|---|---|
| margin `(s₁−s₂)/s₁` (a score, not a probability) | 1.1312 | 0.3974 | 0.4879 | 0.227 |
| softmax, T = 0.03 | 0.4986 | 0.1525 | 0.0905 | 0.771 |
| **platt** | **0.4141** | **0.1358** | **0.0141** | 0.710 |

Platt is adopted as `M`'s probability output. O1's blocker is cleared.

## 2. Arm B at matched pool size

Arm B is subsampled to **exactly 24,517 users**, Arm A's pool size, seed 20261102. Pool size drives
top-1 accuracy and the shape of the score distribution, so comparing a 100,000-user pool with a
24,517-user one would confound pool size with the shock. This is now part of premise P6.

| evaluation users | Arm A (shock) | Arm B (no shock) |
|---|---|---|
| top-1 accuracy | 0.7149 | **0.9111** |
| platt ECE | 0.0141 | **0.0050** |
| platt Brier | 0.1358 | **0.0513** |
| chosen T | 0.03 | 0.03 |

## 3. The drift–score association exists without a shock

| | Arm A (shock) | Arm B (no shock) |
|---|---|---|
| mean JS drift (held-out half) | 0.702 | 0.593 |
| `r`(JS drift, `p_platt`) | −0.657 | **−0.604** |
| `r`(JS drift, correct) | −0.537 | −0.315 |

**This is the result that forces the interpretation to change.** Drift predicts the matcher's
confidence almost as strongly in the control as in the shock arm. The paper may **not** claim that
measured drift predicts linking difficulty only under disruption, nor present the association as a
disruption phenomenon. It is a general property of longitudinal location data.

(The weaker `r`(drift, correct) in Arm B, −0.315 against −0.537, is largely a ceiling effect: Arm B
is 91% accurate, so there are few errors for drift to sort. Arm C's primary configuration gives
−0.423, between the two: [P8_armC.md](P8_armC.md) §4.)

## 4. What the shock does change: calibration

Binned on **absolute** JS drift, not per-arm quantiles. (The per-arm composite used in
`results/P2_decoupling.md` is rank-based, so its quintiles are not comparable between arms — a
methodological error caught and corrected here.) `gap` = accuracy − mean predicted probability, so a
negative gap is overconfidence.

| JS drift bin | A n | A acc | A mean p | A gap | A ECE | B n | B acc | B mean p | B gap | B ECE |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.00–0.45 | 598 | 0.980 | 0.953 | +0.027 | 0.028 | 2187 | 0.984 | 0.986 | −0.002 | 0.003 |
| 0.45–0.55 | 1331 | 0.965 | 0.920 | +0.044 | 0.045 | 2640 | 0.980 | 0.970 | +0.010 | 0.011 |
| 0.55–0.65 | 2334 | 0.919 | 0.853 | +0.066 | 0.066 | 3005 | 0.957 | 0.942 | +0.015 | 0.015 |
| 0.65–0.75 | 3020 | 0.820 | 0.747 | +0.074 | 0.074 | 2524 | 0.885 | 0.881 | +0.005 | 0.014 |
| 0.75–0.85 | 3001 | 0.613 | 0.583 | +0.030 | 0.060 | 1553 | 0.739 | 0.737 | +0.001 | 0.017 |
| **0.85–1.00** | 1975 | 0.218 | 0.460 | **−0.242** | **0.242** | 350 | 0.494 | 0.568 | **−0.074** | 0.089 |

Bootstrap 95% CIs on the top bin (2,000 resamples of evaluation users):

| | n | gap | ECE |
|---|---|---|---|
| Arm A, own calibrator | 1,975 | **−0.242** [−0.261, −0.223] | 0.242 [0.223, 0.261] |
| Arm A, **Arm B's** calibrator | 1,975 | **−0.386** [−0.405, −0.366] | 0.386 [0.367, 0.406] |
| Arm B, own calibrator | 350 | −0.074 [−0.124, −0.027] | 0.089 [0.076, 0.154] |

Three findings.

**(a) Aggregate ECE hides the failure.** Arm A's aggregate ECE is 0.0141, but its
drift-conditional ECE — the coverage-weighted mean of the within-drift-bin ECEs — is **0.0906**, an
understatement of **6.4×**. Arm B understates by 2.8× (0.0050 → 0.0142). The reason is visible in
the table: Arm A's own calibrator is *under*confident at low and middle drift (+0.03 to +0.07) and
*over*confident at the top (−0.24), and the two cancel in the aggregate. **Reporting a single ECE
for an attack under drift is not sufficient**, which is itself a contribution.

**(b) Transport is where the failure is severe.** The realistic deployment case is a calibrator
fitted before a disruption and used during it — nobody can fit on an emergency they have not seen.
Applying Arm B's calibrator to Arm A gives a top-bin gap of **−0.386**: it claims 0.60 probability
where the true rate is 0.22. Its CI does not overlap Arm A's own-calibrator CI. Note that with Arm
B's calibrator the low-drift underconfidence disappears (gaps −0.003 to −0.012), confirming that
Arm A's own-calibrator underconfidence is an artefact of fitting on a shock-contaminated
population.

**(c) At equal measured drift, the shock arm is still harder.** Arm A's accuracy is below Arm B's in
every bin, and the gap widens with drift (0.980 vs 0.984 at the bottom, 0.218 vs 0.494 at the top).
So **measured drift does not fully capture what the disruption does.** This cuts against claim (ii)
being the whole story, and must be stated. It also cannot be explained by Arm A having less data:
Dataset 2 users are 10–13% *more* active than Dataset 1 users
([P4_P10_precheck.md](P4_P10_precheck.md) §2), which should help matching, not hurt it.

## 5. Consequence: the claim in A_PREMISES §0 narrows

Replace "measured drift predicts linking failure under disruption" with:

> Measured drift predicts re-identification difficulty in ordinary conditions as well as under
> disruption. A drift-blind confidence calibrator absorbs ordinary drift (every drift stratum in the
> control is calibrated to within 0.015), but fails in the high-drift stratum a disruption
> populates, and fails far more severely when it is transported from ordinary conditions into the
> disruption. Aggregate calibration metrics understate this by about 6×.

This is narrower than the original phrasing and further from `berrang2026revisiting`, which is good.
It also makes the drift-swap (O3) carry more weight, because the causal question is now about
*calibration transport*, not about whether drift hurts matching.

**Superseded in part, 2026-11-02 (later the same day).** The gap sweep in
[P1_P3_P6_P7.md](P1_P3_P6_P7.md) supplies the mechanism behind §3 and §4, and it should be stated
instead of the correlational framing above: measured drift adds **+0.027 to +0.030 AUC** beyond the
confidence and coverage at *every* gap in Arm A, but **−0.009 to +0.008** — that is, nothing — in
Arm B. In ordinary conditions the score already absorbs drift, so a drift-blind calibrator loses
nothing by ignoring it. Under the disruption drift carries information the score does not, and that
excess is exactly what the calibrator is missing. **Measured drift is diagnostic precisely when the
attack's own score stops being a sufficient statistic for its correctness.**

## 6. What this does not establish

- **One matcher, one profile family.** Nothing shows the calibration failure is matcher-independent.
- **Arms A and B are different user populations**, not the same people with and without a shock.
  Every cross-arm statement is between-population. P4 established they share city, grid and calendar
  period, which is the strongest available control but is not randomisation.
- **No Arm C.** GeoLife (gradual, endogenous drift, raw GPS) is untouched.
- **No causal claim.** §4(b) is an observational transport comparison. O3's drift swap
  ([../B_SWAP_DESIGN.md](../B_SWAP_DESIGN.md)) is what would support a causal reading, and it has
  only been designed and validated, not run.
- The Platt calibrator includes `log(coverage)` as a feature, so the failure is **not** explained by
  coverage: the calibrator already had it.
- The top drift bin in Arm B holds only 350 evaluation users, so its −0.074 gap is the least precise
  number in §4.
- P1 (identity stability) remains unverifiable from the YJMob100K release.
