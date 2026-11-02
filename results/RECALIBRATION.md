# Drift-aware recalibration: the falsifiable prediction, run

Run 2026-11-02 (`../scripts/recalibrate.py`). Discharges the prediction stated in
`../../paper/OUTLINE.md` §7.5 and carried into `main.tex` §VI-D, and the open item in
[`ATYP_BASELINE.md`](ATYP_BASELINE.md) §7 about whether `yuksekgonul2023beyond`'s published
recalibrator repairs what its covariate failed to explain.

> drift-aware recalibration should improve calibration in the shock arm and buy nothing in the
> control. If it helps in the control too, C2b is wrong.

**C2b is not contradicted.** But the pre-registered primary endpoint did **not** come out
cleanly, and the model that worked was a secondary one. Both facts are recorded here in that
order, and §4 explains why the result is weaker than the headline numbers look.

## 1. The pre-registered primary endpoint, reported first

Primary, fixed in the script header before execution: the reduction in **|gap| in the worst
drift quintile** on evaluation users, **AAR(drift) against blind**, per arm, bootstrap 95% CI.
AAR(drift) was chosen as primary because it is `yuksekgonul2023beyond`'s published functional
form with drift substituted for atypicality — a drift-dependent temperature,
`logit' = φ(δ)·logit_blind + b`, `φ(δ) = c₂δ² + c₁δ + c₀`.

| | Arm A (shock) | Arm B (control) |
|---|---|---|
| blind \|gap\| in Q5 | 0.2018 | 0.0103 |
| **AAR(drift) reduction** | **+0.0064 [+0.0016, +0.0111]** | **+0.0094 [−0.0096, +0.0103]** |

The pre-registered rule for reading (a) required Arm A's CI to exclude zero **and** Arm B's
reduction to be less than half of Arm A's point estimate. The first condition holds. **The
second fails**: Arm B's point reduction (+0.0094) is larger than Arm A's (+0.0064).

**So on the pre-registered endpoint the answer is reading (c), the method fails** — not because
drift is uninformative, but because this functional form repairs almost nothing in either arm.
Arm A's 0.2018 gap is still 0.1954 afterwards. Reading (b), C2b contradicted, is **not**
triggered: Arm B's interval contains zero and its baseline gap is 0.0103, so there was nothing
there to repair. See §4 on why the control cannot carry this test.

## 2. The secondary model that does work, labelled as secondary

Adding drift as a **linear feature** alongside the blind logit, rather than as a temperature:

| Arm A (shock), evaluation users | blind | +drift (linear) |
|---|---|---|
| gap by drift quintile | +0.044, +0.069, +0.073, +0.041, **−0.202** | −0.015, −0.010, +0.022, +0.069, **−0.044** |
| drift-conditional ECE | 0.0913 | **0.0397** |
| Brier | 0.1358 | 0.1209 |
| NLL | 0.4141 | 0.3713 |
| AUC | 0.8578 | 0.8885 |
| aggregate ECE | 0.0141 | 0.0157 |
| **Q5 \|gap\| reduction** | — | **+0.1576 [+0.1544, +0.1607]** |

In the control arm the same model gives **−0.0054 [−0.0263, +0.0243]**: the interval contains
zero and the point estimate is slightly *negative*.

That is the predicted asymmetry, and it is large: the worst-quintile gap falls from −0.202 to
−0.044 under the shock and does not improve at all in the control. **But it was not the
pre-registered endpoint**, and §4 says why it does not mean what it appears to.

## 3. Their published method, on their own covariate, recovers nothing

AAR(atyp) is `yuksekgonul2023beyond`'s recalibrator in its own form on its own kind of
covariate — the three rarity estimators from [`ATYP_BASELINE.md`](ATYP_BASELINE.md):

| Arm A | Q5 gap | Q5 \|gap\| reduction | drift-conditional ECE |
|---|---|---|---|
| blind | −0.202 | — | 0.0913 |
| **AAR(atyp)** | **−0.200** | **+0.0022 [+0.0013, +0.0031]** | **0.0917** |
| AAR(drift) | −0.195 | +0.0064 | 0.0868 |
| +drift (linear) | −0.044 | +0.1576 | 0.0397 |

The published method leaves the failure essentially untouched: 0.2018 → 0.1996. This closes the
open item in `ATYP_BASELINE.md` §7 — rarity neither explains the failure nor repairs it — and it
is the strongest single argument in the paper that the covariate had to be a change measure.
AAR(both) adds nothing over AAR(drift), so the two modulators are not complementary.

## 4. Why the headline number is weaker than it looks

**The model that worked is not a recalibrator.** A recalibrator re-maps confidence without
changing who is ranked above whom. `+drift (linear)` raises AUC from 0.8578 to 0.8885, so it
**changes the ranking**: it is a better attack model, not a repair of an existing one. Most of
its apparent calibration gain is C2b restated — drift carries information the score lacks — and
not evidence that a miscalibrated confidence can be fixed after the fact.

The pure rank-preserving recalibrators, AAR(drift) and AAR(atyp), are the honest test of
repairability, and they recover **+0.0064** and **+0.0022** of a 0.2018 gap. The defensible
sentence is therefore:

> Under disruption the matcher's confidence is not merely mis-scaled; the ordering it induces
> is itself wrong on high-drift users, so no drift-dependent temperature can repair it. Drift
> must enter the model, not the calibration layer.

That is a more interesting claim than "drift-aware recalibration works", and it is better
supported. It also explains the aggregate-versus-conditional pattern once more in the repair
direction: AAR(drift) *worsens* Arm A's aggregate ECE from 0.0141 to 0.0425 while improving the
drift-conditional figure, because it distorts the well-calibrated bulk to chase the tail.
Aggregate ECE is the wrong objective here.

**The control cannot carry the falsification as it was written.** Arm B's blind worst-quintile
gap is already 0.0103. A prediction of the form "it should buy nothing in the control" is close
to untestable when the control has nothing to buy: any method scores ≈0 there whether or not it
works. The prediction in `OUTLINE.md` §7.5 was written by me and is weaker than it reads; the
informative half of it is the shock arm alone. A stronger falsification would compare the
*shape* of the repair across arms at matched baseline gaps, which this run does not do.

## 5. What this does not establish

- **No new matcher or drift measure was run.** This only re-maps confidence on the stored
  record, reusing the matcher, the drift values and the calibration/evaluation split.
- **Not tested out of arm.** The practically important question — does a drift-aware model
  fitted on an undisrupted population transfer to a disrupted one? — is the transported-calibrator
  experiment, not this one. This fits and evaluates inside the same arm.
- **One functional form each.** `φ` is quadratic in the modulator because that is the published
  choice; an isotonic or spline recalibrator conditioned on drift was not tried and might do
  better than +0.0064.
- **Arm C is untouched**, at n = 56 calibration is not assessable there.

## 6. Reproduce

```
python3 scripts/recalibrate.py results/cal_armA_dataset2.npz results/atyp_armA_dataset2.npz armA_dataset2
python3 scripts/recalibrate.py results/cal_armB_dataset1.npz results/atyp_armB_dataset1.npz armB_dataset1
```

Console transcript: `recalibrate_output.txt`. Artefacts: `recal_arm{A,B}_*.npz`, excluded from
git because they carry per-user probabilities.
