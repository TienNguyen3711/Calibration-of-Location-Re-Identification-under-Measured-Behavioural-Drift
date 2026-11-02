# Resolving premise P2: drift / score decoupling

Run 2026-11-02, Arm A only (YJMob100K Dataset 2). Settles the premise that
[../A_PREMISES.md](../A_PREMISES.md) §4 calls "the central threat":

> **P2** If `δ_u` is computed from the same statistics that `M` compares, then a large `δ_u`
> *implies* a poor match by construction, and "drift predicts failure" is circular rather than
> empirical.

**Verdict: P2 is discharged, with one residual caveat recorded in §6.** The association between
measured drift and the matcher's confidence survives on records the matcher never sees, it is
person-specific, and the circularity it would have carried is quantified rather than assumed.


> **Correction, 2026-11-02.** The first pass computed Spearman with an argsort-of-argsort
> implementation that ignores tied ranks. The outcome `y` is binary and the `anchor` component is
> binary, so every correlation involving either was attenuated — by 0.10 to 0.19 for `r(·, y)`. All
> such numbers below are the corrected, midrank values (`scipy.stats.spearmanr`). Correlations
> between two continuous variables were unaffected, so every ECE, calibration gap, ΔAUC and
> bootstrap interval in these documents is unchanged. `scripts/matchlib.py`, `p2_decoupling.py`,
> `calibrate.py` and `compare_arms.py` were all fixed and the affected analyses re-run.

Full console output: [P2_decoupling_output.txt](P2_decoupling_output.txt). Variables:
`p2_vars_dataset2.npz`.

Reproduce (from `Codebase/code/`):
```
gzcat data/raw/yjmob100k/yjmob100k-dataset2.csv.gz | awk -f scripts/p2_profiles.awk > prof_d2.tsv
python3 scripts/p2_decoupling.py prof_d2.tsv results/p2_vars_dataset2.npz
```
`prof_d2.tsv` is a 99,126-line intermediate (~100 MB) and is not kept in the repo; the awk pass
takes about 36 s and the analysis about 8 s. Harmless `RuntimeWarning: ... in matmul` lines come
from this machine's Accelerate BLAS; every saved variable was checked finite.

## 1. The chosen scheme: interleaved-day split

P2's obligation allowed two routes — disjoint features, or a held-out temporal split within each
window. **The held-out split is now fixed as the design**, in this specific form:

Within every window, days are assigned by parity of the day index:

| Half | Days | Consumed by |
|---|---|---|
| **M** | even day index | the matcher `M` (profiles `π_u`) |
| **D** | odd day index | the drift measure `δ_u` |

No record contributes to both halves. `W_ref = d 45..59` gives 7 M-days and 8 D-days; `W_tgt =
d 60..74` gives 8 M-days and 7 D-days.

**Why parity and not first-half / second-half.** A contiguous split would make the drift measurement
cover a different calendar sub-period from the matcher, confounding drift with sub-period. Parity
keeps both halves spanning the full window. It also balances weekday almost exactly: over a 15-day
window the even-day set hits each of the seven weekday phases once (`d mod 7` = 4,6,1,3,5,0,2 for
d = 46..58) and the odd-day set hits each once with one repeat. Weekday mix is therefore not a
confound between the two halves.

**What the split does and does not remove.** It removes *record sharing* — the mechanical
construction P2 names. It does not remove *behavioural dependence*: a person's adjacent days are
autocorrelated, so the two halves describe the same person in the same fortnight. That dependence is
the phenomenon under study, not leakage.

## 2. Matcher as instrumentation

Profiles are 500 m cell-visit histograms over the M half, square-root (Hellinger) transformed,
IDF-weighted, L2-normalised; score = cosine; pool = every user with all four (window × half) groups.

| | |
|---|---|
| users in pool | **24,517** of 24,999 |
| top-1 accuracy | **0.7128** |
| top-10 accuracy | **0.9230** |
| median rank of the true user | 1 |

This is well inside the range where confidence is worth calibrating: a matcher near 0 or near 1
would leave nothing to study. Per A_PREMISES §0 the accuracy is a reported precondition, not a
contribution.

## 3. The decoupling result

Spearman correlations, n = 24,517. `c` = margin confidence `(s₁−s₂)/s₁`, `y` = top-1 correct.
**Decoupled** = δ from the D half (the design). **Coupled** = the same δ formula computed on the
matcher's own M half (the circular counterfactual, computed only to be reported).

| drift component | r(δ_dec, c) | r(δ_cou, c) | r(δ_dec, y) | r(δ_cou, y) | r(dec, cou) | r(δ_dec, coverage) |
|---|---|---|---|---|---|---|
| visit-distribution JS | **−0.415** | −0.475 | **−0.537** | −0.602 | 0.917 | −0.218 |
| cell-set dissimilarity | −0.339 | −0.406 | −0.434 | −0.509 | 0.765 | −0.125 |
| anchor change | −0.135 | −0.151 | −0.222 | −0.239 | 0.524 | −0.033 |
| novelty rate | −0.344 | −0.375 | −0.436 | −0.475 | 0.864 | −0.230 |
| entropy change | −0.086 | −0.100 | −0.128 | −0.139 | 0.351 | −0.041 |
| slot-histogram JS | −0.120 | −0.123 | −0.163 | −0.166 | 0.719 | −0.484 |
| **composite** | **−0.391** | −0.443 | −0.508 | −0.563 | 0.845 | −0.303 |

Three readings.

**(a) No leakage.** P2 warned that "a near-deterministic relationship is evidence of leakage, not of
a finding". The strongest decoupled correlation is **−0.415**. Substantial, nowhere near
deterministic.

**(b) The circularity premium is small but real, and now measured.** Coupled drift correlates with
`c` at −0.475 against the decoupled −0.415. In the nested fit (§4) the inflation is clearer: the
apparent AUC gain from adding drift is **+0.062** decoupled and **+0.077** coupled, so roughly
**20%** of what a naive same-records analysis would have claimed is record sharing alone. Reporting
this number is the substance of discharging P2.

**(c) The split costs almost no signal.** `r(dec, cou) = 0.917` for the visit-distribution
component. This is a test–retest reliability estimate from two disjoint day sets: drift is a stable
person-level property, so holding out half the days weakens the measure very little. It also frames
(a) properly — a measure with reliability 0.92 correlating with `c` at only −0.42 is clearly not the
same quantity as `c`.

## 4. Preliminary nested fit (this is O2, shown only as a smell test)

Logistic fit for `y`; `c`, `log(coverage)` and drift all z-scored.

| model | mean log-lik | AUC | Brier |
|---|---|---|---|
| `c` alone | −0.4476 | 0.8289 | 0.1480 |
| `c` + coverage | −0.4426 | 0.8335 | 0.1463 |
| `c` + coverage + drift **decoupled** | **−0.3632** | **0.8952** | **0.1186** |
| `c` + coverage + drift **coupled** | −0.3383 | 0.9103 | 0.1103 |

Drift adds beyond confidence and coverage on held-out records. This is **not** the paper's O2 result
— no held-out users, no regularisation, no calibration-specific metric, Arm A only, one matcher —
and must not be cited as one. It is here because a null at this stage would have meant the design
had nothing to measure.

**Placebo.** Giving each user another user's drift, with their own score kept, collapses the
correlation to zero: 20 permutations give mean `r(δ, c)` = **+0.0008** (sd 0.0059, max |r| 0.0128)
and mean `r(δ, y)` ≈ 0 to the same precision. The association is specific to the person, which is what the O3
drift-swap will test properly.

**The margin score is not a probability.** Mean `c` = 0.227 against accuracy 0.713, so its raw
Brier (0.396) is not a miscalibration measure and is not reported as one. The deciles are monotone
and strongly informative (accuracy 0.377 in the lowest decile, 0.946 in the 0.3–0.4 decile), so the
score carries the ordering a calibrator needs. **O1 requires `M` to emit a probability**; adding
that mapping is the next build step, not a finding.

## 5. Consequences for the design

**C4 — Fix the primary drift measure as visit-distribution JS on the held-out half.** It is the
strongest (−0.415), by far the most reliable (0.917), and only mildly coverage-confounded (−0.218).

**C5 — Demote slot-histogram JS.** `r(δ, coverage) = −0.484`, stronger than its correlation with the
score (−0.120). It mostly measures how much data a person has. Keep it as a reported component with
coverage adjustment; never in the composite without one.

**C6 — Anchor change and entropy change are weak here** (−0.135, −0.086), and anchor change has
the lowest reliability of the spatial components (0.524). This agrees with C1 from the P4/P10
pre-check. Both stay as reported components, neither drives the composite.

**C7 — Coverage must be a covariate, not a control variable only.** `r(coverage, c) = +0.168` and
`r(composite drift, coverage) = −0.303`: coverage predicts both, so it enters every model (P7).

## 6. What this does not establish, and the residual caveat

- **Residual caveat on P2.** The decoupling is **record-level, not feature-level**. The strongest
  drift component (visit-distribution JS) is the same *family* of statistic as the matcher's
  profile, computed on disjoint records. The one genuinely feature-disjoint component
  (slot-histogram JS, temporal where the matcher is spatial) correlates with `c` at only −0.120 and
  is coverage-confounded. So the headline association is strongest where drift and profile share a
  feature family. This is not circularity in P2's sense — no shared records — but it must be stated
  in the paper, and the feature-disjoint component must be reported beside the primary one so a
  reader can see the difference. Risk **R3** in `gaps.md` moves from *open* to *mitigated and
  documented*, not to *closed*.
- Arm A only. The same decoupling must be rerun on Arm B (Dataset 1) and Arm C (GeoLife) before any
  claim is made; Arm B in particular must show whether the drift–score association exists without a
  shock.
- One matcher. Nothing here shows the result is matcher-independent.
- No calibration result. `M` does not yet emit a probability, so O1 (calibration loss under drift)
  has not been touched.
- P1 (identity stability) remains unverifiable from the YJMob100K release.
