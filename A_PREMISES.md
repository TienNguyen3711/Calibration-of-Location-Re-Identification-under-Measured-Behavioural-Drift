# A: premises for the drift-calibration experiment

Status: design contract, written 2026-11-01 **before any experiment code exists**. It fixes what
will be measured and what would invalidate the measurement. It is **not** a result, and it does not
claim any premise below is already satisfied. Every premise carries an obligation and a current
status; most statuses are "not yet implemented".

Data provenance, licences and the audit of the raw files are in
[data/README.md](data/README.md) and [data/LICENCES.md](data/LICENCES.md). Positioning and the gap
this work targets are in `../../Material/_ai_research/intake/gaps.md`.

## 0. The claim, and what is deliberately not claimed

**Claim (narrowed 2026-11-02 on measured evidence; see
[results/O1_calibration.md](results/O1_calibration.md) §5).** On location data, under per-person
behavioural drift that is *measured* rather than assumed, (i) a drift-blind confidence calibrator
absorbs ordinary drift but loses calibration in the high-drift stratum a disruption populates, and
far more severely when it is transported from ordinary conditions into the disruption; (ii) measured
drift predicts which confident matches are wrong, over and above the elapsed time gap, the attack's
own confidence and observability; (iii) aggregate calibration metrics understate this by about 6×,
so calibration under drift must be reported drift-conditionally; and (iv) it is the drift
actually present in the target window, rather than a fixed property of the individual, that makes a
transported confidence calibrator overconfident — tested by a drift-swap counterfactual.

*Claim (iv) reworded 2026-11-02 on measured evidence.* The earlier wording, "the relationship is
specific to each person's own drift", is **contradicted** by the amended O3 analysis: a stranger's
drift does the same damage as the person's own, a random donor is as damaging as a nearest
neighbour, and once applied drift is known the recipient's own drift carries no information
(b = −0.0001 [−0.0071, +0.0075]). The pre-registered O3 verdict is nonetheless **INCONCLUSIVE**, and
the reworded claim rests on an exploratory amendment with one control still missing
([results/O3_swap.md](results/O3_swap.md) §7).

**Mechanism (measured 2026-11-02).** Drift is diagnostic precisely when the attack's own score
stops being a sufficient statistic for its correctness: measured drift adds ΔAUC +0.026…+0.030
beyond confidence and coverage at every gap in the shock arm, and ≈0 in the control. That is why a
drift-blind calibrator stays calibrated throughout the control and fails in the shock arm's
high-drift stratum.

*Superseded phrasing, kept for the record:* "the confidence score of a re-identification attack
loses calibration" and the implication that the drift–score association is a disruption phenomenon.
Measured 2026-11-02: the association is nearly as strong in the no-shock control (`r` = −0.604
against −0.657), so it is a general property of longitudinal location data, not a shock effect.

**Not claimed.**

- No claim that our matcher is more accurate than MainTUL, Lu 2026, Tournier 2022 or any other.
  The matcher is instrumentation. Its accuracy is a precondition to be reported, not a contribution.
- No claim that "drift degrades re-identification" is new. That is established in another modality
  (`berrang2026revisiting`). The contribution is confined to calibration, on location data, with a
  swap test.
- No claim of field-wide novelty. The gap statement is **about the 32-source corpus only** until the
  R1 search in `gaps.md` is done.
- No causal claim about *why* a person drifts. Drift is a covariate here, not an outcome.
- **No claim that the drift–score association is specific to disruptions.** It is not
  (§0 above). Only the *calibration* failure is disruption-specific.
- No claim that measured drift explains everything the disruption does. At equal measured drift the
  shock arm is still less linkable in every stratum, so the drift measure is incomplete
  ([results/O1_calibration.md](results/O1_calibration.md) §4c).

## 1. Experimental arms and their distinct roles

Three arms over two independent sources. The arms are **not** three replications.

| Arm | Data | Role | Question it alone answers |
|---|---|---|---|
| **A** | YJMob100K Dataset 2 (25,000 users; d 0..74, shock at d 60..74) | Primary | Does measured drift predict confident-but-wrong matches at a fixed gap? |
| **B** | YJMob100K Dataset 1 (100,000 users; d 0..74, no shock) | **Control** | Is the Arm-A effect attributable to the shock, or to window position / passage of time alone? |
| **C** | GeoLife 1.3 (~40-80 usable users, multi-year) | **Replication** under a different drift regime | Does the finding survive gradual, endogenous drift and raw GPS? |

Only **A and C** are independent evidence. **B is a control**: it shares Arm A's city, resolution,
collection process and window geometry, and differs in the one factor of interest. Reporting B as a
third replication would misstate the design.

Measured 2026-11-01 (`results/P4_P10_precheck.md`): B is a **stronger** control than this paragraph
originally assumed. Dataset 1 and Dataset 2 cover not only the same city and grid but the same
calendar period and weekday phase — their daily mobility series correlate at **+0.990** across
days 0..59, and that correlation falls to **−0.336** across days 60..74. Arm B therefore holds city,
grid, season, weekday structure and window position fixed at once.

Arms never share individuals: see P1. No cross-dataset linking occurs anywhere in this project.

## 2. Window geometry (fixed in advance; do not retune after seeing outcomes)

Let `W_ref` and `W_tgt` be disjoint day intervals, `|W_ref| = |W_tgt| = L`, and `g` the number of
days between the end of `W_ref` and the start of `W_tgt`.

**Arms A and B — primary configuration.** `L = 15`, `W_tgt = d 60..74`, `W_ref = d 45..59`, `g = 0`.

The gap is zero by construction, so in Arm A **drift varies across people while the gap is held
constant**. That is the cleanest possible form of the "beyond the gap" test, because gap cannot
explain any between-person variation.

**Arms A and B — secondary gap sweep.** Same `L`, same `W_tgt`, `W_ref` moved earlier:
`d 30..44` (g = 15), `d 15..29` (g = 30), `d 0..14` (g = 45).

**Arm C.** `L` in days, chosen so that enough users retain coverage in both windows; `g` swept.
The usable user count must be computed by the build script, not assumed (see P7).

## 3. Definitions

- **Profile** `π_u(W)`: the representation of user `u` over window `W` that the matcher consumes.
- **Matcher** `M`: given a target profile and a candidate pool `P`, returns a ranking and a
  confidence score `c_u ∈ [0,1]` for its top candidate. Fixed, known design (profile/histogram
  family, after `tournier2022expanding`); replaceable without changing the rest of the pipeline.
- **Outcome** `y_u ∈ {0,1}`: 1 if the top-1 candidate is the true user.
- **Drift** `δ_u`: a scalar (and a small vector of components) measuring how much user `u`'s
  behaviour changed between `W_ref` and `W_tgt`.

**Drift is not a single quantity.** Candidate components, each already motivated by a card in the
corpus, to be computed per user per window pair and reported separately as well as combined:

| Component | Definition sketch | Source of the idea |
|---|---|---|
| visit-distribution distance | L1 / Jensen-Shannon between cell-visit distributions | `dong2021statistical` |
| anchor change | does the inferred home / most-visited cell change | `golle2009anonymity`, `stroebl2025investigating` |
| top-N set overlap | Jaccard of the N most-visited cells | `qi2019user` |
| stationarity / regularity change | difference in the two scalars between windows | `teixeira2021impact` |
| entropy change | difference in per-window location entropy | `teixeira2021impact` |
| novelty rate | share of target-window cells never seen in reference | `mcinerney2013breaking` |

`δ_u` must be **computed at a fixed spatial and temporal resolution**, because these quantities move
strongly with resolution (`teixeira2021impact`).

**Record split (fixed 2026-11-02, P2).** `π_u` is built from the **even** day indices of each window
and `δ_u` from the **odd** ones. Every component above is computed on the odd half only. The same
formulas on the even half are computed as well, but solely as the circularity counterfactual that
P2 obliges us to report; they never enter a result.

**Measured drift properties (2026-11-02, Arm A, n = 24,517).** Visit-distribution JS is the primary
scalar: it is the strongest predictor of the matcher's confidence (Spearman −0.415), the most
reliable (test–retest 0.917 across the two day halves), and only mildly confounded by coverage
(−0.218). Slot-histogram JS is **demoted**: its correlation with coverage (−0.484) exceeds its
correlation with the score (−0.120), so it is reported only with coverage adjustment and is excluded
from the composite until adjusted. Anchor change (−0.111, reliability 0.563) and entropy change
(−0.086) carry little signal here, consistent with the noise floors above; both remain reported
components and neither drives the composite.

**Measured noise floors (2026-11-01, Arm B, no shock).** Between two adjacent 15-day windows with no
disruption at all, **30.1%** of users already change their night-time home cell and **24.2%** change
their top cell, and the mean Jaccard of their visited-cell sets is only **0.258**. Under the shock
these move to 50.0%, 50.7% and 0.198. Two consequences:

- *Anchor change* carries real signal (about +20 points) but has a baseline near a coin flip at
  500 m resolution, where a genuine move and an estimate flipping to an adjacent cell are
  indistinguishable for one user. It must be reported next to its control baseline and must not be
  the primary drift measure.
- *Set-overlap* measures start from a low base (about 28 of 69 reference cells revisited), so
  profiles built on set intersection will be weak here irrespective of drift. This supports the
  distribution/histogram profile family already chosen for `M`.

## 4. Premises, with obligations

**P1 — Identity stability.** Within one dataset, a `uid` denotes the same person for the whole
observation period. If ids were recycled or re-pseudonymised mid-period, "drift" would partly be
identity change, and the whole measurement collapses.
*Obligation:* state the assumption explicitly; it is **not verifiable** from the YJMob100K release,
which carries no documentation on this point. Cross-dataset `uid` identity is separately **known to
fail** (Dataset 1 and Dataset 2 both number from zero: see `data/raw/yjmob100k/SOURCE.md`), so no
analysis may link them.

**P2 — Drift and score must not be mechanically coupled (the central threat).** If `δ_u` is computed
from the same statistics that `M` compares, then a large `δ_u` *implies* a poor match by
construction, and "drift predicts failure" is circular rather than empirical.
*Obligation:* `δ_u` is computed on a feature set **disjoint** from the matcher's profile, or on a
**held-out temporal split within each window** that the matcher never sees. The chosen scheme must
be stated, and the empirical correlation between `δ_u` and `c_u` reported. A near-deterministic
relationship is evidence of leakage, not of a finding. This is risk R3 in `gaps.md`.

*Scheme chosen and implemented, 2026-11-02 — **interleaved-day split**.* Within every window, days
with an **even** index feed the matcher (`π_u`) and days with an **odd** index feed `δ_u`; no record
enters both. Parity rather than a contiguous split, so both halves span the whole window and the
weekday mix is balanced (each of the 7 weekday phases appears once in the even-day set of a 15-day
window). Measured: `r(δ, c) = −0.415` on held-out days against −0.475 on the matcher's own days, so
the circularity the split removes is about 20% of the apparent predictive gain; drift test–retest
reliability across the two halves is 0.917; the placebo with another user's drift gives
`r ≈ +0.001`. Residual caveat: the split is record-level, not feature-level — see
[results/P2_decoupling.md](results/P2_decoupling.md) §6. R3 is **mitigated and documented**, not
closed.

**P3 — Drift must vary at fixed gap.** The claim is "beyond the gap". If `δ` and `g` are collinear,
the claim is untestable.
*Obligation:* the primary configuration fixes `g = 0`, so between-person variation in `δ` cannot be
gap. For the gap sweep, report `corr(δ, g)` and fit the nested models within gap strata.

*Discharged 2026-11-02, both arms.* `r(δ, g)` = +0.016 (A) and +0.074 (B); within-gap sd of drift
exceeds the between-gap sd of the means by 31× and 11×. The nested fits within gap strata produced
the paper's mechanism: drift adds **ΔAUC +0.026…+0.030 at every gap in Arm A and ≈0 in Arm B**, so
in ordinary conditions the score already absorbs drift, and only under disruption does drift carry
information the score lacks. **Caveat:** in Arm A the target window is always the shock window, so
that arm's gap sweep is confounded and the gap reading must be taken from Arm B.
[results/P1_P3_P6_P7.md](results/P1_P3_P6_P7.md)

**P4 — Control comparability (Arm B).** Arm B is a valid control only if Dataset 1 and Dataset 2
share city, grid, resolution and broadly similar population behaviour, differing in the shock alone.
*Obligation:* **unverified.** Compare cell-occupancy distributions and basic mobility statistics
between the two files before using B as a control. If they differ materially, B is a weak control
and must be reported as such.

**P5 — Swap plausibility.** Transplanting user `v`'s drift onto user `u` must yield a counterfactual
that is still a plausible mobility profile. Location profiles carry geographic and continuity
constraints that feature vectors do not (`berrang2026revisiting` swapped feature vectors, which are
unconstrained).
*Obligation:* define the swap on a representation where the operation is closed — for example,
applying `v`'s *change* in visit distribution relative to `v`'s own anchors to `u`'s anchors, rather
than copying raw cells. Validate swapped profiles before interpreting results: if the swap produces
implausible traces, "swap restores accuracy" is an artefact of implausibility, not evidence of
person-specificity. Three conditions to compare, after `berrang2026revisiting`: own drift,
nearest-neighbour drift, population-predicted drift.

**P6 — Pool definition.** Confidence depends on the candidate pool. A score of 0.9 against 100
candidates and against 25,000 are not the same quantity.
*Obligation:* fix pool size and composition in advance, identically across drift strata, and report
sensitivity to pool size. Any pool that varies with `δ` invalidates the calibration comparison.

**P7 — Observability versus behaviour.** A user can appear to drift because they were *observed*
less, not because they *moved* differently. This is severe for GeoLife (median user: 19 active days)
and real for YJMob100K (presence in a 30-minute slot is itself a behaviour).
*Obligation:* compute per-user coverage in each window; filter on a pre-declared coverage rule;
include coverage as a covariate in the nested models; report how many users each rule removes.
Drift measures that are sensitive to sample size (entropy in particular) must be estimated on
equalised sample sizes or bias-corrected.

*Discharged 2026-11-02, with one obligation we could not meet.* The coverage curve is reported in
full rather than a threshold being chosen, because outcomes had already been inspected by the time
the rule would have been declared — a chosen threshold would be post-hoc selection dressed as
pre-registration. The primary analysis therefore applies **no coverage filter**, with
`log(coverage)` as a calibrator feature. The result is flat against the filter (Arm A top-bin gap
−0.257 → −0.269 when 41% of users are dropped). Entropy change is **not** the sample-size trap this
premise assumed (`r` with coverage −0.040 in A, −0.083 in B); `js_slots` is (−0.484).

**P8 — Representation for Arm C.** To run one protocol on both sources, GeoLife must be discretised
onto a comparable grid. Cell size and time bin change every drift quantity.
*Obligation:* fix the discretisation a priori (target: ~500 m cells, 30-minute slots, to match
YJMob100K), justify it, and report a sensitivity check at one coarser and one finer setting.

*Discharged 2026-11-02.* 500 m / 30-minute slots, one observation per user per slot, timestamps
shifted GMT→UTC+8 (the `.plt` files are GMT; without the shift the daily rhythm is displaced by
eight hours). Sensitivity at 250 m and 1000 m: `r`(drift, correct) = −0.378, −0.423, −0.480, every
CI excluding zero, while top-1 moves 0.750 → 0.732 → 0.607. **The deeper finding is a
representation mismatch the premise did not anticipate:** GeoLife logs *trips*, YJMob100K samples
*presence*, so the same discretisation gives a GeoLife user a median of only **6.7 occupied bins per
active day out of 48**. Usable users were computed, not assumed: 56 (relative anchoring, L = 56,
g = 0, ≥40 bins) against 36 (calendar-aligned, L = 90). **Arm C cannot speak to calibration** at
that pool size — see P9 and [results/P8_armC.md](results/P8_armC.md) §5.

**P9 — Estimation precision.** Reliability curves and ECE need enough observations per bin. Arm C
has tens of users.
*Obligation:* pre-specify binning; report bootstrap CIs over users; do not state calibration
conclusions for Arm C that the sample size cannot support. Arm C is corroboration, not a precise
estimate.

**P10 — The shock must be real.** Arm A assumes d 60..74 of Dataset 2 is a genuine behavioural
disruption. Its nature is undisclosed, and some users may be unaffected.
*Obligation:* verify the shift empirically at population level before contrasting A with B. The
per-person design does not depend on every user being affected — that heterogeneity is the point —
but if the shock is invisible in aggregate, Arm A has no contrast with Arm B.

**P11 — Ethics and terms.** YJMob100K is CC BY 4.0 with no re-identification clause; GeoLife names
location privacy as an intended use; the HuMob 2024 four-city data **forbid re-identification** and
are quarantined.
*Obligation:* no script reads `data/raw/_QUARANTINE_humob2024_4cities/`. Re-identification is
pseudonym-to-pseudonym within a release and cannot reach a natural person; state this in the paper's
ethics section together with the YJMob100K authors' "prediction benchmark" intended-use statement.
See `data/LICENCES.md`.

## 5. Outputs (fix before running; do not add outcomes afterwards)

**Precondition, reported but not claimed.** Matcher top-1 / top-10 accuracy per arm, so a reader can
judge whether the attack is strong enough for conclusions about it to matter.

**O1 — Calibration under drift.** Reliability curve and ECE of `c_u` against `y_u`, computed overall
and **stratified by drift decile**. Brier score per stratum. False-discovery rate among
high-confidence matches (`c_u > 0.95` and `> 0.9`), per stratum — the quantity that matters
operationally, following `rocher2019estimating` and `tournier2022expanding`.

**O2 — Incremental value of drift.** Nested logistic regression for `y_u`:
baseline `c_u + g + coverage` versus `+ δ_u`. Report the likelihood-ratio test, ΔAUC and ΔBrier,
with bootstrap CIs over users. This is the direct test of the "beyond gap and confidence" claim.

**O3 — Swap counterfactual.** O1 and O2 recomputed under own drift, nearest-neighbour drift and
population-predicted drift, with the plausibility validation from P5.

**O4 — Component breakdown.** Which drift components in §3 carry the signal, reported honestly
including components that carry none.

**Reporting rule added 2026-11-02.** Every calibration number must be reported **drift-conditionally
as well as in aggregate**, with the ratio between the two stated. An aggregate ECE alone is not an
acceptable report: in Arm A it understates the drift-conditional ECE by 6.4× because underconfidence
at low drift cancels overconfidence at high drift. O1 must also include the **transport** case — a
calibrator fitted on the control and applied to the shock arm — which is the realistic deployment
scenario and where the failure is worst (top-bin gap −0.386).

## 6. Audit: no premise silently assumed satisfied

| Premise | How it will be checked | Status 2026-11-02 |
|---|---|---|
| P1 identity stability | Not checkable in general; the recycling failure mode is now bounded by a break-pattern test against a shuffled-pair null | **Assumption, now bounded** 2026-11-02 — break pattern in 0.075% of control users (1.21% in the shock arm, attributable to the shock: 66.8% of its weakest links sit at the shock boundary against 26.3% in the control). [results/P1_P3_P6_P7.md](results/P1_P3_P6_P7.md) |
| P2 drift/score decoupling | Disjoint features or held-out split; report `corr(δ, c)` | **DISCHARGED** 2026-11-02 (Arm A) — interleaved-day split implemented; `r(δ, c) = −0.415` decoupled vs −0.475 coupled; circularity premium ≈ 20% of the AUC gain; placebo ≈ 0. One residual caveat (record-level, not feature-level decoupling). [results/P2_decoupling.md](results/P2_decoupling.md) |
| P3 drift varies at fixed gap | `g = 0` primary design; `corr(δ, g)` on the sweep | **PASS** 2026-11-02 both arms — `r(δ, g)` = +0.016 (A) and +0.074 (B); within-gap sd of drift exceeds between-gap sd of the means by 31× (A) and 11× (B). Drift adds ΔAUC +0.026…+0.030 in every gap stratum of A and ≈0 in B — the mechanism behind O1 |
| P4 control comparability | Cell-occupancy comparison, Dataset 1 vs Dataset 2 | **PASS** 2026-11-01 — same city *and* same calendar period; daily-series corr +0.990 over d 0..59. [results/P4_P10_precheck.md](results/P4_P10_precheck.md) |
| P5 swap plausibility | Closed-form swap on anchor-relative change + validation | **Run 2026-11-02; O3 INCONCLUSIVE** under the pre-registered rule ([results/O3_swap.md](results/O3_swap.md)). Self-swap exact at full scale (24,517 users); S2 NN clears the gate (50.9% admissible, 0.80 delivery) but the decision test is null on both terms; S1/S2 disagree, which is the declared inconclusive trigger. Also **designed and machinery validated** 2026-11-02 — [B_SWAP_DESIGN.md](B_SWAP_DESIGN.md); self-swap identity check PASS (p95 JS 1e-8, anchor exact 3000/3000); plausibility gate declared; binding limitation found and recorded (≈⅓ of swaps inadmissible, 80% magnitude delivery). No O3 result |
| P6 pool definition | Fixed pool size, sensitivity sweep | **DISCHARGED** 2026-11-02 — arms compared at an identical pool; sweep over 1k/2.5k/5k/10k/23.8k in both arms. The failure **grows** with pool size (top-bin gap −0.137 → −0.252 in A), so the headline is conservative. Pool 1,000 in Arm B is degenerate (98.6% accuracy, calibrator unidentified) and must not be quoted; minimum usable pool ≈ 2,500 |
| P7 observability | Coverage per window, pre-declared filter, covariate | **DISCHARGED** 2026-11-02 — per-user coverage computed; the filter is reported as a curve, **not chosen** (its pre-declaration obligation could not be met; see P7 below). Discarding 41% of Arm A users moves the top-bin gap only −0.257 → −0.269. Entropy change turns out **not** coverage-driven (−0.040); `js_slots` is (−0.484). Earlier partial result stands: during the shock Dataset 2 has +4.5% more users observed per day but −9.8% records. GeoLife counts in `SOURCE.md` |
| P8 Arm C representation | Fixed grid + sensitivity at two other settings | **DISCHARGED** 2026-11-02 — 500 m / 30 min fixed a priori, GMT→UTC+8 corrected, sensitivity at 250 m and 1000 m. Usable users **computed**: 56 (relative, L=56, g=0, ≥40 bins) and 36 (calendar, L=90). `r`(drift, correct) = −0.423 [−0.620, −0.193] primary, negative with CI excluding zero in 6 of 7 settings. **Calibration not assessable at n=56.** [results/P8_armC.md](results/P8_armC.md) |
| P9 estimation precision | Pre-specified bins, bootstrap CIs | **In use** 2026-11-02 — absolute drift bins declared in `scripts/compare_arms.py`; bootstrap 95% CIs on the top-bin gap reported ([results/O1_calibration.md](results/O1_calibration.md) §4) |
| P10 shock is real | Population-level shift check on Dataset 2 | **PASS** 2026-11-01 — between-file corr falls +0.990 → −0.336 at d 60; cells/user-day −12.4% against +1.0% in the control |
| P11 ethics and terms | Quarantine honoured; ethics section written | **Drafted** 2026-11-02 — [../paper/ETHICS.md](../paper/ETHICS.md); quarantine enforced in the repo (no script reads that path), not only on paper. Three open author actions: GeoLife citations, courtesy notice, institutional-review question |

## 7. Acceptance criteria before the claim may be stated

1. ~~P2 is implemented and the drift/score correlation is reported.~~ **Satisfied for Arm A
   2026-11-02** — see [results/P2_decoupling.md](results/P2_decoupling.md). Still to do before the
   claim: rerun the same decoupling on Arm B and Arm C, and report the feature-disjoint drift
   component beside the primary one (§6 caveat there).
2. ~~P4 and P10 are checked.~~ **Satisfied 2026-11-01** — both pass; see
   [results/P4_P10_precheck.md](results/P4_P10_precheck.md). Arm B's role is confirmed and may be
   stated more strongly than first drafted.
3. P5's swap is defined, implemented and validated for plausibility. An unvalidated swap cannot
   support the causal part of the claim. **Satisfied as to the method, 2026-11-02; the result is
   INCONCLUSIVE.** Defined and pre-registered in [B_SWAP_DESIGN.md](B_SWAP_DESIGN.md), machinery
   validated by an exact self-swap at full scale, plausibility gate declared with measured pass
   rates, donor ladder run. **Consequence: the causal part of §0 claim (iv) may not be stated.**
   See [results/O3_swap.md](results/O3_swap.md) §3 for why the pre-registered analysis plan could
   not detect the effect, and §4 for the declared amendment.
4. P7's coverage rule is declared before outcomes are inspected, and the number of users it removes
   is reported for every arm.
5. The R1 search in `gaps.md` is done, so the novelty statement can move from "in this corpus" to a
   defensible claim — or be narrowed.
6. Arm C is reported with its sample size attached to every statement. **Satisfied 2026-11-02** —
   every Arm C number in [results/P8_armC.md](results/P8_armC.md) carries its `n`, and §5 there
   states plainly that Arm C corroborates the drift–outcome association only and says nothing about
   the calibration findings.

Until all six hold, results are internal and the claim in §0 stays conditional.

## 8. Change log

- 2026-11-01 — created, before any experiment code. Window geometry, drift components, premises
  P1-P11 and outputs O1-O4 fixed. No code written yet.
- 2026-11-01 — **P4 and P10 pre-checked and both pass** (`scripts/precheck_p4_p10.awk`,
  `results/P4_P10_precheck.md`; every row of both files read). Three design consequences folded back
  in: Arm B is a stronger control than drafted (§1); anchor-change drift has a 24-30% noise floor and
  cell-set overlap a low base (§3); P7 is partly discharged for Arm A, since the shock shows more
  users observed but fewer records. Acceptance criterion 2 satisfied. P2 remains the open risk.
- 2026-11-02 — **P2 discharged for Arm A** (`scripts/p2_profiles.awk`, `scripts/p2_decoupling.py`,
  `results/P2_decoupling.md`). The interleaved-day split is now part of the design (§3, §4-P2); the
  first working matcher exists as instrumentation (cell-visit Hellinger/IDF cosine, top-1 0.713,
  top-10 0.923 on a 24,517-user pool). Four design consequences folded back in: visit-distribution
  JS is the primary drift scalar, slot-histogram JS is demoted for coverage confounding, anchor and
  entropy change are weak, and coverage enters every model. Acceptance criterion 1 satisfied for
  Arm A. Open: the same check on Arms B and C, a probability output for `M` (needed for O1), and
  P5's swap.
- 2026-11-02 — **`M` now emits a probability; Arm B run; P5 designed** (`scripts/calibrate.py`,
  `scripts/compare_arms.py`, `scripts/swap_feasibility.py`, `results/O1_calibration.md`,
  `B_SWAP_DESIGN.md`). Platt on score order statistics plus `log(coverage)` is adopted as `M`'s
  probability output (Arm A eval ECE 0.0141, Brier 0.1358); a top-K softmax shortcut was rejected
  because it misses the exact full-pool value by up to 0.76. **§0 narrowed on measured evidence:**
  the drift–score association is nearly as strong in the no-shock control, so it is not a disruption
  phenomenon; only the calibration failure is. Two further findings folded in: aggregate ECE
  understates drift-conditional ECE by 6.4× in Arm A, and at equal measured drift the shock arm is
  still less linkable, so the drift measure is incomplete. P6 and P9 partly discharged; acceptance
  criterion 3 partly satisfied. Open: Arm C, the O3 donor ladder, and the R1 novelty search.
- 2026-11-02 — **P1, P3, P6, P7 and P11 closed out; the paper's mechanism found**
  (`scripts/windows_profiles.awk`, `scripts/matchlib.py`, `scripts/sweeps.py`,
  `scripts/p1_identity.py`, `results/P1_P3_P6_P7.md`, `../paper/ETHICS.md`, `../paper/OUTLINE.md`).
  All five 15-day windows are extracted, so `g` sweeps over 0/15/30/45 in both arms at a matched
  pool of 23,763. **§0 gains a mechanism**: drift adds ΔAUC +0.026…+0.030 beyond confidence and
  coverage at every gap in Arm A and ≈0 in Arm B, so drift is informative exactly when the score
  stops being sufficient. P1's recycling bound tightened to ≈0.1% using the control. P3 PASS in both
  arms, with the Arm A confound recorded. P6 discharged with the warning that the failure grows with
  pool size and that pool 1,000 is degenerate. P7 discharged with an explicit statement that its
  pre-declaration obligation could not be met. Day 27 exposure quantified (4.4% of records, at
  `g = 30` only). **P8 (Arm C / GeoLife) is now the only untouched premise.**
- 2026-11-02 — **P8 discharged; all eleven premises now resolved** (`scripts/ingest_geolife.py`,
  `scripts/geolife_feasibility.py`, `scripts/armc_run.py`, `results/P8_armC.md`). GeoLife ingested
  to 84,741 occupied 30-minute bins from 24,876,975 point lines; representation fixed a priori at
  500 m / 30 min with GMT→UTC+8 correction and sensitivity at 250 m and 1000 m. Usable users
  computed: 56 primary, 36 in the protocol-identical calendar variant. Arm C **corroborates** the
  drift–outcome association (−0.423 [−0.620, −0.193], and negative with CI excluding zero in 6 of 7
  configurations) but **cannot assess calibration** at that pool size, so the replication claim is
  scoped to C1/C4 only. Two honest negatives recorded: the calendar configuration's CI contains
  zero, and GeoLife has no shock, so Arm C replicates the Arm B regime rather than Arm A's.
- 2026-11-02 — **O3 run; INCONCLUSIVE; and an amendment to B_SWAP_DESIGN §7, declared here as
  that document requires.** The donor ladder was run on all 24,517 Arm A users
  (`scripts/o3_swap.py`, `scripts/o3_analysis.py`, `results/O3_swap.md`). The swap machinery is
  verified: the self-swap reproduces the real world exactly at full scale. S2 with a
  nearest-neighbour donor clears both pre-registered gate thresholds, yet the decision test returns
  **null on both applied and own drift** (b = −0.0026 [−0.0079, +0.0028] and +0.0011
  [−0.0046, +0.0066], VIF 1.0, n = 6,224). S1 cannot separate the two (collinear by construction,
  VIF 2.5), and S1/S2 disagreeing in direction is the declared inconclusive trigger. **A 3,000-user
  pilot had shown a significant effect in the predicted direction; it did not replicate at full
  scale, and the interim report of it was wrong.** Diagnosis: §7's instruction to refit the
  calibrator inside the swapped world absorbs the applied drift — per-bin gaps are ≈0 in every S2
  bin while accuracy falls from 0.767 to 0.264 — so the plan removed the transport phenomenon it was
  meant to test.
  **AMENDMENT (date 2026-11-02, reason above):** `scripts/o3_swap.py` now also transports the real
  world's Platt calibrator into every swapped world unchanged and reports its ECE beside the
  refitted one. The amended analysis is **exploratory and reported separately**; the pre-registered
  verdict stays INCONCLUSIVE regardless of what it shows.
  Three design costs now measured and to be reported: per-half reference anchors agree for only
  **70.2%** of users (a cost of the P2 decoupling, unanticipated); **6.6–6.8%** of synthetic mass
  leaves the 200×200 grid under S2; and the **POP drift-field rung is degenerate** (pooled novelty
  mass 0.007, so applying the population field nearly erases drift and top-1 rises to 0.9991) —
  withdrawn as designed, with the S1 population-magnitude variant as its substitute. That
  degeneracy contrasts with `berrang2026revisiting`, where a population-level drift prediction
  restored attack AUC to 0.821, and the contrast needs explaining in the paper.
- 2026-11-02 — **Amended O3 run; §0 claim (iv) reworded.** Transporting the real world's calibrator
  into the swapped worlds raises S2 ECE from 0.0127 to **0.1609** and makes the decision test
  decisive: the residual tracks **applied** drift (b = −0.0625 [−0.0688, −0.0559] for an NN donor,
  −0.0633 [−0.0700, −0.0571] for a random donor) and **not** the recipient's own drift
  (−0.0001 [−0.0071, +0.0075]), with a monotone per-bin gap of −0.083 → −0.255. This **contradicts
  the old wording of claim (iv)** and supports a stronger causal reading with weaker
  person-specificity; §0 reworded accordingly. Still exploratory: the amendment followed the
  pre-registered result, S1 continues to disagree with S2, and **generic profile disruption is not
  excluded** — S2 cuts top-1 from 0.711 to 0.494, and no perturbation control exists that moves the
  profile as much without changing measured drift. That control must be built before the causal
  claim is printed.
- 2026-11-02 — **Correction to every tie-affected correlation.** The Spearman helper used
  argsort-of-argsort, which ignores tied ranks; the binary outcome `y` and the binary `anchor`
  component are all ties, so those correlations were attenuated by 0.10–0.19. Fixed to midranks in
  `scripts/matchlib.py`, `p2_decoupling.py`, `calibrate.py`, `compare_arms.py` and the affected
  analyses re-run. `r`(JS drift, correct) becomes −0.537 (A) and −0.315 (B), from −0.432 and −0.127, and `r(δ, g)`
  becomes +0.016 (A) and +0.074 (B) because `gap` takes only four values.
  **No ECE, calibration gap, ΔAUC or bootstrap interval changed**, because those never involve a
  binary variable in a rank correlation. Documents corrected in place with a correction note.

- **2026-11-02, the atypicality baseline (new obligation, created and discharged the same day).**
  Reading Yuksekgonul et al., NeurIPS 2023 (`10.52202/075280-1669`) in full showed that **ECE by
  per-sample atypicality quantile is already published** (their Table 1), so the *shape* of O1 §4b
  is not ours, and the clause "our covariate is not derived from the model" is **false** and was
  deleted — their class-atypicality estimator needs only training labels
  (`Material/_ai_research/intake/R1_ATYPICALITY_RULING.md`). That created a mandatory baseline:
  does a plain rarity score of the target window alone explain the calibration failure? **It does
  not** (`results/ATYP_BASELINE.md`). Primary endpoint ΔAUC(drift | base + atypicality) =
  **+0.0560 [+0.0509, +0.0612]** under the shock against **+0.0043 [+0.0002, +0.0087]** in the
  matched control; three rarity estimators spanning their menu each add ≈0 on their own in both
  arms; the overconfidence collapse (gap −0.202, ECE 0.2018 in the worst drift quintile) has no
  counterpart on any rarity axis, whose largest gap is +0.019 and is *under*-confidence. Rarity and
  drift are near-orthogonal with opposite signs: `r(a_pop, y) = +0.082`, `r(js_cells, y) = −0.537`,
  `r(a_pop, js_cells) = +0.003`. Reported rather than hidden: `atypicality | base + drift` is
  +0.0245, a suppression effect between drift and the score-space `a_knn` (r = +0.709), which is a
  real out-of-sample gain but not an exogenous rarity signal. Nothing was refitted; user order was
  verified against the stored run at `max|Δ| = 0` in both arms, and the C2b number reproduces
  (+0.0307 here against +0.0292 published at g = 0).

- **2026-11-02, the O3 control (S4) and the second matcher family.** The missing
  equal-magnitude, drift-preserving control was built and run (`results/S4_CONTROL.md`). The M
  half is perturbed by Dirichlet mixing over each user's own support, bisected to match S2-NN's
  top-1 (0.4966 against 0.4943), while the D half is left real so measured drift is unchanged by
  construction and verified at `max|d_app - d_own| = 0.000e+00`. Transported-calibrator ECE
  **0.0383** against S2-NN's **0.1609**, ratio **0.238**, pre-registered threshold 0.5:
  **SURVIVES**. Generic profile disruption no longer explains the O3 §7 result. Two
  qualifications are recorded rather than smoothed: the sweep's worst point reaches 0.477 of S2's
  figure, so the separation is twofold at worst, not fourfold; and the secondary check I
  pre-registered in that script was **mis-specified** -- I justified it by "`d_app` is constant by
  construction", which is false (`d_app` equals each user's own varying drift), so its failure at
  b = -0.0413 is the paper's main effect reappearing inside the control, not a broken
  construction. The property that rule meant to protect is verified separately and exactly.
  O3's own pre-registered verdict stays INCONCLUSIVE, since S1 still disagrees with S2.
