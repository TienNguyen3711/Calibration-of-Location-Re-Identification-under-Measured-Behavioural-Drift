# B: the drift-swap counterfactual (premise P5, output O3)

Status: **design fixed and machinery validated 2026-11-02. No O3 result exists.** This document
is the pre-registration: it fixes the construction, the donor ladder, the plausibility gate and the
decision rule *before* any swap is run through the matcher, so the conclusion cannot be chosen after
seeing it.

Companion to [A_PREMISES.md](A_PREMISES.md). Validation script:
[scripts/swap_feasibility.py](scripts/swap_feasibility.py); output
[results/P5_swap_feasibility_output.txt](results/P5_swap_feasibility_output.txt).

## 0. What the swap is for

[results/O1_calibration.md](results/O1_calibration.md) §4 shows that a confidence calibrator fitted
in ordinary conditions becomes badly overconfident in the high-drift stratum of a disruption. That
is observational: high-drift people might simply *be* different people, and drift a marker rather
than a cause.

The swap asks the counterfactual directly: **if this person had drifted like someone else, would
their confidence have been miscalibrated to the degree the donor's drift predicts, or to the degree
their own drift predicts?** The permutation placebo already run
([results/P2_decoupling.md](results/P2_decoupling.md) §4) tests whether the *statistic* is
person-specific; it does not construct the counterfactual world, so it cannot answer this.

## 1. Representation

For user `w` and window `W`, let `a_w(W)` be the **anchor**: the modal cell of the M half of `W`.
Behaviour is then described by the distribution `q_w(W)` over **offsets** `o = cell − a_w(W)`.

Drift splits into two transportable parts:

| part | definition |
|---|---|
| anchor displacement | `D_w = a_w(W_tgt) − a_w(W_ref)` |
| offset-distribution change | `q_w(W_ref) → q_w(W_tgt)` |

The anchor is the **modal** cell, not the night-time home cell, because at 500 m the night-home
estimate flips to an adjacent cell for 30% of users between adjacent undisturbed windows
([results/P4_P10_precheck.md](results/P4_P10_precheck.md) §4, C1) and is therefore too unstable to
be a reference frame.

**Measured (Arm A, 3,000 users).** The frame is only moderately concentrated: median mean offset
radius **11.0 cells**, with a median of **0.43** of the reference mass within 3 cells of the anchor
and **0.63** within 10. The frame is usable but it is not a tight local frame, so any claim that the
swap preserves "the person's neighbourhood" is limited to that 63%.

## 2. The two swaps

### S2 — transport swap (primary)

`u` receives donor `v`'s drift:

```
anchor    ã_u(tgt) = a_u(ref) + D_v
offsets   w(o) = q_u,ref(o) · min(C, q_v,tgt(o) / q_v,ref(o))        for o in supp(q_v,ref)
          w(o) = q_u,ref(o) · κ                                      otherwise, κ = mean ratio
          ~q = (1 − m_nov) · w / Σw   +   m_nov · (novel offsets)
novelty   m_nov = Σ q_v,tgt(o) over o outside supp(q_v,ref)
coverage  total record count held at u's own real target count
```

`C = 50` caps the ratio. **No additive smoothing is used**, deliberately: with `ε > 0` the
self-swap identity check (§3) fails by construction, and an identity check that cannot pass is
worthless.

**Novelty placement — two modes.** The donor's novel offsets cannot simply be reused, or a donor
would relocate `u` into the donor's own geography, destroying the thing the swap is supposed to hold
fixed. So:

- `novelty="exact"` reuses the donor's own novel offsets. **Used only for the §3 identity check.**
- `novelty="radius"` resamples the same number of novel offsets at matched radii, in uniformly random
  directions, on offsets `u` has not visited. **This is the analysis mode.** Direction information
  is discarded; that is a stated loss.

### S1 — magnitude-only swap (secondary)

`~q = (1−λ)·q_u,ref + λ·q_u,tgt`, with `λ` solved by bisection so that `JS(q_u,ref, ~q)` equals the
JS that donor `v` actually showed. This carries the donor's drift **size** but `u`'s own drift
**direction**. It is far more plausible by construction — it never leaves the geodesic between two
real distributions of the same person — but it is a weaker counterfactual. It is the fallback if S2
fails the plausibility gate.

## 3. Correctness: the self-swap identity check

With `v = u` the support condition `supp(q_u,ref) = supp(q_v,ref)` holds, so the ratio is exact and
S2 must reproduce `u`'s real target distribution exactly. **Tolerance: p95 `JS(synthetic, real) <
1e-6` and the anchor reproduced exactly.**

**Result: PASS.** Median JS **2.9 × 10⁻⁹**, p95 **1.0 × 10⁻⁸**, anchor exact for 3,000/3,000.
**4 users of 3,000 (0.13%)** exceed 1e-6, in every case because the ratio cap `C = 50` binds where a
donor offset has tiny reference mass and large target mass. Those users are reported, never silently
dropped.

In analysis mode (`novelty="radius"`) the self-swap reproduces the real target only to median
JS ≈ 0.65, by construction — novel locations are resampled. So §3 validates the **ratio and anchor
machinery**; it does not validate novelty placement, which only the plausibility gate constrains.

## 4. The donor ladder

Mapping A_PREMISES §5 O3 ("own drift, nearest-neighbour drift, population-predicted drift") onto
concrete donors, weakest counterfactual first:

| rung | donor | what it tests |
|---|---|---|
| **S0** | drift *label* permuted, no synthesis | is the statistic person-specific (done: `r ≈ +0.001`) |
| **S3** | `v = u` | machinery correctness (done: PASS) |
| **NN** | nearest neighbour on baseline cell count, reference entropy and `log(coverage)` | a similar person with different drift — the conservative donor, and the primary |
| **RND** | uniformly random donor | upper bound on implausibility |
| **POP** | the population-mean drift field in place of any donor | does any individual drift signature matter, or only the average |

## 5. The plausibility gate (P5's obligation)

A synthetic profile is **admissible** only if all four marginals fall inside the central 95% band of
the **real** target-window profiles:

| marginal | real 95% band (Arm A) |
|---|---|
| distinct cells | [11, 99] |
| offset entropy (bits) | [1.880, 6.285] |
| largest-cell share | [0.061, 0.620] |
| mean offset radius (cells) | [1.083, 48.197] |

*Implementation requirement:* the band above is **unconditional** and is the preliminary, looser
gate. The final gate must be **matched on coverage decile**, because all four marginals move with
how much data a person has.

**Measured pass rates and transfer fidelity, 3,000 users:**

| swap | donor | all four marginals | `r`(applied drift, donor's drift) | magnitude delivered |
|---|---|---|---|---|
| S1 magnitude | random | 68.0% | 0.625 | 89% |
| S1 magnitude | **nearest neighbour** | 67.8% | **0.821** | **92%** |
| S2 transport | random | 65.2% | **0.908** | 80% |
| S2 transport | **nearest neighbour** | 63.8% | 0.899 | 80% |

## 6. The obstacle found, and what it costs

**Drift at this resolution is mostly novelty, and novelty is what makes a synthetic profile
implausible.** A donor places a median **0.56** of its target mass on offsets it never visited in
its reference window, while only a median **0.29** of `u`'s reference mass sits on offsets the donor
also has. So the multiplicative ratio field acts on under a third of `u`'s mass, and over half the
transported drift arrives as novelty.

Controlling the synthetic cell count to keep it inside the real band was **tried and rejected**: it
necessarily discards donor novelty mass, and transfer fidelity collapsed from `r = 0.92` to
`r = 0.38`. **The two goals are in direct conflict** — you cannot both deliver the donor's drift
magnitude and hold the cell count inside the real band, because the drift *is* the new cells. This
is recorded as the binding limitation of P5 rather than engineered away.

Consequences, all of which must appear in the paper:

1. **About a third of swaps are inadmissible.** O3 is reported on the admissible subset only.
2. **That subset is selected on drift**, since high-drift swaps fail the gate more often. The O3
   analysis must therefore be stratified on applied drift, or weighted back to the full
   distribution, and the drift distribution of the dropped set reported beside the kept one.
3. **S2 delivers only 80% of the donor's drift magnitude.** Any O3 effect size is an underestimate
   by roughly that factor, and must be presented as a lower bound.
4. **Direction information is discarded** by radius-matched novelty placement, so S2 transports
   drift *magnitude and radial shape*, not drift *direction*. A donor-direction-preserving variant
   is possible but would reintroduce the geography-relocation problem; it is not in scope.

## 7. Decision rule, fixed before running

Run the ladder, recompute `M`'s probability (`scripts/calibrate.py`'s Platt model, refitted on
calibration users of the swapped world) and recompute the drift-stratified calibration table of
[results/O1_calibration.md](results/O1_calibration.md) §4. Then:

- **Claim (iii) is supported** if, on admissible swaps, the drift-stratified calibration gap tracks
  the **applied** drift and not the recipient's own drift — specifically, if regressing the gap on
  applied drift and own drift jointly leaves applied drift significant and own drift not.
- **Claim (iii) fails** if the gap tracks the recipient's own drift, or their identity, irrespective
  of the applied drift. Drift would then be a *marker* of hard-to-link people, not a cause of
  calibration failure, and the paper must say so and drop the causal part of §0.
- **O3 is inconclusive, not negative**, if the admissible fraction falls below 50%, if S2's
  magnitude delivery falls below 60%, or if S1 and S2 disagree in direction. An inconclusive O3 is
  reported as inconclusive; it is not resolved by relaxing the gate.

Nothing in this document may be changed after the first swap is run through the matcher. Changes
go in A_PREMISES §8 with a date and a reason.

**AMENDED 2026-11-02, after the first run.** The instruction above to recompute `M`'s probability
by refitting the Platt model "on calibration users of the swapped world" is wrong: refitting inside
the counterfactual lets the calibrator learn the swapped world's difficulty, absorbing the applied
drift so the decision test has nothing to measure (per-bin gaps ≈ 0 in every S2 bin while accuracy
falls 0.767 → 0.264). The amended run additionally **transports the real world's calibrator
unchanged** into every swapped world, mirroring [results/O1_calibration.md](results/O1_calibration.md)
§4b where the effect lives. The amendment is recorded in `A_PREMISES.md` §8 with its reason; the
pre-registered verdict for the first run stands as **INCONCLUSIVE**
([results/O3_swap.md](results/O3_swap.md)), and the amended analysis is exploratory.
