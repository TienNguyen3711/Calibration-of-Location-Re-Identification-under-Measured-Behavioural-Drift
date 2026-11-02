# O3: the drift-swap counterfactual — result

Run 2026-11-02 on Arm A (YJMob100K Dataset 2), all 24,517 users with all four
(window × half) groups. Executed against the rule fixed in
[../B_SWAP_DESIGN.md](../B_SWAP_DESIGN.md) §7 **before** any swap was run.

Scripts: [../scripts/o3_swap.py](../scripts/o3_swap.py),
[../scripts/o3_analysis.py](../scripts/o3_analysis.py). Output:
[O3_swap_prereg_output.txt](O3_swap_prereg_output.txt) (the swap run),
[O3_swap_output.txt](O3_swap_output.txt) (the decision test). Variables: `o3_vars.npz`,
`o3_result.pkl`.

## Verdict: **INCONCLUSIVE** under the pre-registered rule, and the reason is diagnosable

The rule said O3 is inconclusive if the admissible fraction falls below 50%, if S2's magnitude
delivery falls below 60%, **or if S1 and S2 disagree in direction**. The third trigger fires, and
S2-RND also misses the admissibility floor by 0.8 points.

Claim (iv) — that the relationship is specific to each person's own drift — is **not supported and
not refuted**. The decision test returned null on both terms in the only condition able to
separate them.

## 1. The ladder, at full scale

| condition | n | admissible | top-1 | ECE (eval) | mean `d_app` | `d_app`/`d_donor` | `r(d_app, d_own)` |
|---|---|---|---|---|---|---|---|
| REAL (no swap) | 24,517 | 88.9% | 0.7113 | 0.0086 | 0.753 | — | +1.000 |
| **S3 self, exact** | 24,517 | **88.9%** | **0.7113** | 0.0090 | **0.753** | — | +1.000 |
| S2 NN donor | 24,441 | **50.9%** | 0.4943 | 0.0127 | 0.602 | **0.80** | +0.334 |
| S2 RND donor | 24,431 | **49.2%** | 0.4931 | 0.0132 | 0.599 | 0.80 | −0.026 |
| S2 POP field | 24,517 | **0.0%** | 0.9991 | 0.0005 | 0.118 | — | +0.450 |
| S1 NN donor | 24,517 | 67.5% | 0.7311 | 0.0071 | 0.693 | 0.92 | +0.791 |
| S1 POP magnitude | 24,517 | 58.4% | 0.7427 | 0.0118 | 0.710 | — | +0.847 |

**The machinery is correct.** S3 (donor = the user, exact novelty) reproduces the real world on
every quantity: top-1 0.7113 against 0.7113, mean `d_app` 0.753 against 0.753, admissible 88.9%
against 88.9%, and a per-bin calibration table identical to within one user per bin. The residual
ECE difference (0.0090 against 0.0086) comes from synthetic counts being floats where the originals
are integers, which perturbs `sqrt(v/tot)` in the last bits and flips a handful of users. At full
scale this is the strongest possible evidence that the swap is implemented as specified.

**Gate outcomes.** S2 NN clears both thresholds (50.9% ≥ 50%, delivery 0.80 ≥ 0.60). S2 RND misses
admissibility at 49.2%. **The anchors are estimated separately per half to keep the halves disjoint,
and they agree for only 70.2% of users** — a cost of the P2 decoupling that was not anticipated and
must be reported.

**The POP rung is degenerate and must be withdrawn as designed.** The pooled population offset
distribution has support almost everywhere, so its novelty mass is **0.007**: applying the
population drift *field* to an individual nearly **erases** drift (mean `d_app` 0.118 against 0.753
real) and pushes top-1 to 0.9991. Nothing passes the plausibility gate. This is itself a finding —
**an individual's drift bears almost no resemblance to the population-average drift field** — and it
contrasts with `berrang2026revisiting` N30, where a population-level drift *prediction* restored
attack AUC to 0.821, close to the nearest-neighbour 0.804. The likely explanation is that location
data carries geographic structure that an unconstrained feature vector does not. The S1 POP
*magnitude* variant is the usable substitute.

## 2. The decision test

OLS of the calibration residual `y − p` on z-scored applied and own drift, admissible evaluation
users only, bootstrap 95% CI over 2,000 resamples.

| condition | n | b(`d_app`) | 95% CI | b(`d_own`) | 95% CI | VIF |
|---|---|---|---|---|---|---|
| **S2 NN donor** | 6,224 | **−0.0026** | [−0.0079, +0.0028] | **+0.0011** | [−0.0046, +0.0066] | **1.0** |
| S2 RND donor | 6,022 | −0.0027 | [−0.0075, +0.0023] | −0.0017 | [−0.0074, +0.0039] | 1.0 |
| S1 NN donor | 8,270 | +0.0433 | [+0.0310, +0.0560] | −0.0609 | [−0.0719, −0.0504] | 2.5 |
| S1 POP magnitude | 7,190 | +0.1088 | [+0.0680, +0.1494] | −0.1351 | [−0.1775, −0.0924] | 33.0 |

**S2 is the only condition that can answer the question, and it answers "neither".** Applied and own
drift are nearly independent there (`r` = 0.334, VIF 1.0), so the design worked — and both
coefficients are null with tight intervals.

**S1 cannot be read as support.** Its `d_app` is the recipient's own drift rescaled to the donor's
magnitude, so the two regressors are collinear by construction (`r` = 0.791, VIF 2.5). Two
significant coefficients of opposite sign under collinearity is the signature of **one effect split
between two correlated predictors**, not of applied drift mattering positively. S1 POP magnitude,
at VIF 33.0, is uninterpretable and is reported only for completeness.

**A smoke-test result that did not replicate, recorded because it was reported.** On a 3,000-user
pilot, S2 NN gave b(`d_app`) = −0.0328 [−0.0521, −0.0133], significant and in the direction claim
(iv) predicts, with b(`d_own`) null. At 24,517 users (n = 6,224 admissible) it is −0.0026
[−0.0079, +0.0028]. **The pilot finding was noise and the interim report of it was wrong.**

## 3. Why the test was null: the pre-registered analysis plan destroyed the effect

The per-bin tables in [O3_swap_output.txt](O3_swap_output.txt) §3 show it directly. In S2 NN the
calibration gap is approximately zero in **every** absolute-drift bin:

| `d_app` bin | n | accuracy | mean p | gap |
|---|---|---|---|---|
| 0.00–0.45 | 2,074 | 0.767 | 0.771 | −0.004 |
| 0.45–0.55 | 1,491 | 0.642 | 0.653 | −0.011 |
| 0.55–0.65 | 1,391 | 0.451 | 0.455 | −0.004 |
| 0.65–0.75 | 933 | 0.318 | 0.347 | −0.028 |
| 0.75–0.85 | 326 | 0.264 | 0.255 | +0.009 |

Accuracy falls from 0.767 to 0.264 across the bins and **the predicted probability tracks it almost
exactly**. In the real world the same table is monotone and signed (+0.011, +0.034, +0.056, +0.044,
−0.009, −0.039).

B_SWAP_DESIGN §7 instructed that `M`'s probability be "refitted on calibration users of the swapped
world". **That instruction is the error.** Refitting inside the counterfactual lets the calibrator
learn the swapped world's difficulty, so the applied drift is absorbed and nothing is left for the
decision test to measure. The phenomenon O1 found is a **transport** failure — a calibrator from one
regime applied to another ([O1_calibration.md](O1_calibration.md) §4b). By refitting, O3 removed
exactly what it was built to test.

This is a flaw in the plan I wrote, not in the swap, and not a property of the data.

## 4. Amendment, declared rather than slipped in

B_SWAP_DESIGN §7 states that nothing in it may change after the first swap is run, and that changes
go in `A_PREMISES.md` §8 with a date and a reason. Accordingly:

> **Amendment, 2026-11-02.** `scripts/o3_swap.py` now additionally carries the **real world's**
> Platt calibrator into every swapped world unchanged, and reports its ECE beside the refitted one.
> Reason: the pre-registered refit absorbs the applied drift (§3 above).

The transported variant is reported **separately** and is **not** the pre-registered result. The
pre-registered verdict stands as INCONCLUSIVE whatever the amendment shows; the amended analysis is
exploratory and must be labelled so in the paper.

## 5. What this means for the paper

- **C5 is a method contribution with a documented negative**, which is what
  `Material/_ai_research/intake/R1_NOVELTY_VERDICT.md` anticipated. The swap is implemented,
  validated by an exact self-swap at full scale, and pre-registered; its first application was
  inconclusive for a reason now identified.
- **The causal part of §0 claim (iv) cannot be stated.** Until the amended analysis is read, the
  paper claims only the observational transport result from O1, and says plainly that the
  counterfactual test was inconclusive.
- **Three design costs are now measured and must be reported**: per-half anchors agree for only
  70.2% of users; 6.6–6.8% of synthetic mass leaves the 200×200 grid in S2; and roughly half of S2
  swaps are inadmissible, with the kept set selected on drift.
- **S2 degrades matching heavily** (top-1 0.711 → 0.494). The swapped world is not merely a
  re-drifted version of the real one; it is a harder problem. Any future iteration must report this
  rather than treat the two worlds as comparable.

## 6. What this does not establish

- Nothing about whether drift *causes* calibration failure. That question remains open.
- The amended, transported analysis was not available when this verdict was fixed.
- Arm A only; no Arm B or Arm C swap was run.
- The plausibility gate is the unconditional-marginal test matched on coverage decile; a stricter
  gate would shrink the admissible set further.

---

# 7. The amended analysis (exploratory): the transported calibrator

Run 2026-11-02 after the amendment in §4. **This is not the pre-registered result**, which stays
INCONCLUSIVE. It is labelled exploratory everywhere it appears.

## 7.1 Transporting the calibrator restores the phenomenon

| condition | ECE, refitted in the swapped world | ECE, real-world calibrator **transported** |
|---|---|---|
| REAL (no swap) | 0.0086 | 0.0086 *(identical by construction)* |
| S3 self, exact | 0.0090 | 0.0088 *(sanity check passes)* |
| **S2 NN donor** | 0.0127 | **0.1609** |
| **S2 RND donor** | 0.0132 | **0.1574** |
| S1 NN donor | 0.0071 | 0.0117 |
| S1 POP magnitude | 0.0118 | 0.0175 |

Refitting inside the counterfactual gave ECE 0.0127; carrying the real world's calibrator into it
unchanged gives **0.1609 — 12.7× worse**. §3's diagnosis is confirmed: the pre-registered plan was
measuring a calibrator that had already learned the swapped world.

## 7.2 The decision test, on the transported probability

Same OLS of `y − p` on z-scored applied and own drift, admissible evaluation users, 2,000 bootstrap
resamples. `d_app` and `d_own` are nearly independent in S2 (`r` = 0.334, VIF 1.0).

| condition | n | b(`d_app`) | 95% CI | b(`d_own`) | 95% CI |
|---|---|---|---|---|---|
| **S2 NN donor** | 6,224 | **−0.0625** | **[−0.0688, −0.0559]** | **−0.0001** | **[−0.0071, +0.0075]** |
| **S2 RND donor** | 6,022 | **−0.0633** | **[−0.0700, −0.0571]** | −0.0051 | [−0.0119, +0.0018] |
| S1 NN donor | 8,270 | +0.0436 | [+0.0308, +0.0562] | −0.0621 | [−0.0732, −0.0504] |
| S1 POP magnitude | 7,190 | +0.1142 | [+0.0734, +0.1572] | −0.1367 | [−0.1808, −0.0946] |

And the gap is monotone in applied drift, reproducing O1 §4b inside the counterfactual:

| `d_app` bin, S2 NN | n | accuracy | mean p | gap | ECE |
|---|---|---|---|---|---|
| 0.00–0.45 | 2,074 | 0.767 | 0.850 | −0.083 | 0.083 |
| 0.45–0.55 | 1,491 | 0.642 | 0.779 | −0.137 | 0.137 |
| 0.55–0.65 | 1,391 | 0.451 | 0.648 | −0.197 | 0.197 |
| 0.65–0.75 | 933 | 0.318 | 0.573 | −0.255 | 0.255 |
| 0.75–0.85 | 326 | 0.264 | 0.505 | −0.241 | 0.241 |

**In both S2 conditions the transported calibrator's miscalibration tracks the drift that was
applied, and the recipient's own drift adds nothing — its coefficient is −0.0001.** Two independent
donor schemes give coefficients that agree to the third decimal (−0.0625, −0.0633).

## 7.3 This supports a causal reading of drift, and **contradicts claim (iv) as §0 words it**

A_PREMISES §0 (iv) says the relationship is *"specific to each person's own drift"*. The amended
analysis says the opposite of that literal wording:

- transplanting **a stranger's** drift produces the same calibration failure as the person's own;
- a **random** donor is as damaging as a **nearest-neighbour** donor (−0.0633 against −0.0625);
- once the applied drift is known, **the person's own drift carries no information** about whether
  the transported calibrator is wrong about them.

The defensible reading is therefore **stronger on causality and weaker on person-specificity**:

> It is the drift actually present in the target window that makes a transported confidence
> calibrator overconfident, not a fixed property of the individual. Drift is not merely a marker of
> hard-to-link people.

That is the useful claim, and it is what `B_SWAP_DESIGN.md` §7's rule calls "supported" — the rule
was written as "the gap tracks the applied drift and not the recipient's own drift". **§0's wording
and §7's rule were never equivalent, and the data has now separated them.** §0 must be reworded;
§7's operationalisation is the one that was pre-registered and the one to keep.

**Contrast with `berrang2026revisiting` (N30), stated carefully.** There, substituting the
nearest neighbour's drift *restored* attack AUC from 0.557 to 0.804, which the authors read as the
degradation being tied to each person's own drift. Our outcome is different — calibration residual,
not attack power — so the two are not directly comparable, and the contrast must be presented as
"different outcome variable, possibly different phenomenon", not as a contradiction of their result.

## 7.4 Why this is still exploratory, and the alternative explanation that is not excluded

- **It is an amendment made after seeing the pre-registered result.** No amount of internal
  consistency changes that. The paper reports the pre-registered verdict first.
- **S1 still disagrees with S2** (+0.0436 against −0.0625). S1's regressors are collinear by
  construction (VIF 2.5), so its signs are not interpretable, but the declared inconclusive trigger
  "S1 and S2 disagree in direction" **still fires even under the amendment**.
- **Difficulty is confounded with drift.** S2 cuts top-1 from 0.7113 to 0.4943, so the transported
  calibrator meets a world that is harder, not only differently drifted. The residual `y − p_tr` is
  a function of the matcher's scores, which absorbs some of this, but **the test does not control
  for the matcher's score features directly.** A variant that adds them as controls is the obvious
  next step and was not run.
- The nearly identical NN and RND coefficients can be read two ways: as evidence that *applied*
  drift is what matters, or as evidence that **any** large perturbation of the profile breaks the
  transported calibrator, drift-shaped or not. Distinguishing these needs a perturbation that
  changes the profile as much as a swap does **without** changing measured drift — a control that
  does not currently exist in this design and should be built before the causal claim is made in
  print.

Until that control exists, the honest statement is: **the amended analysis is consistent with
applied drift causing transported-calibrator failure, and rules out the recipient's own drift as the
explanation, but does not yet rule out generic profile disruption.**

---

## 8. Generic profile disruption is now ruled out (added 2026-11-02)

The control this section called for exists: [`S4_CONTROL.md`](S4_CONTROL.md). Perturbing the
matcher's input until top-1 falls to **0.4966**, matching S2-NN's 0.4943, while leaving the
drift half untouched (`max|d_app - d_own| = 0.000e+00`), costs a transported-calibrator ECE of
**0.0383** against S2-NN's **0.1609** -- a ratio of **0.238** against a pre-registered threshold
of 0.5. **Verdict: SURVIVES.** The collapse is not explained by the swapped world merely being
harder.

Two qualifications carried from that file: at the worst point of the lambda ladder a
structureless perturbation reaches 0.0767, which is **0.477** of S2's figure, so the separation
is about twofold rather than fourfold; and the secondary check pre-registered in `s4_control.py`
failed because I specified it on a false premise, which S4_CONTROL.md §4.2 states in full.

The pre-registered verdict of **this** section is unchanged and remains **INCONCLUSIVE**: it was
triggered by S1 and S2 disagreeing in direction, which S4 does not address.
