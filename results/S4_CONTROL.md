# S4: the equal-magnitude, drift-preserving control

Run 2026-11-02 (`../scripts/s4_control.py`). Closes the hole that
[`O3_swap.md`](O3_swap.md) §7.4 names as the largest in the study, and that
`../../paper/OUTLINE.md` listed as a blocking experiment.

**Pre-registered verdict: SURVIVES.** An equally damaging, equally large perturbation of the
matcher's input, carrying **no drift structure at all**, does not reproduce the transported
calibrator's failure: ECE **0.0383** against S2-NN's **0.1609** at matched top-1 accuracy, a
ratio of **0.238** against the pre-registered threshold of 0.5.

Two things qualify that, both in §4: the margin is a factor of ~2 rather than ~4 if one takes the
worst point of the λ sweep instead of the accuracy-matched point, and **the pre-registered
secondary check failed — because I specified it on a false premise**, not because the
construction is broken.

## 1. What was missing, and what this is

O3 §7.4: the transport swap cuts top-1 from 0.711 to 0.494, so the swapped world is not only
differently drifted but **harder**, and a random donor damages the transported calibrator as much
as a nearest neighbour. So perhaps *any* large perturbation breaks a transported calibrator,
drift-shaped or not. Excluding that needs a perturbation that moves the profile as much as a swap
without changing measured drift.

**Construction.** Drift is measured on the odd-day (D) half; the matcher consumes the even-day
(M) half. S4 perturbs the **M half only** and leaves the D half exactly as in the real world, so
measured drift is unchanged by construction. The perturbation mixes the target M-half offset
distribution toward a random Dirichlet direction **over that user's own support**:

    mix = (1 − λ)·q_target + λ·Dirichlet(1,…,1)

The anchor is not displaced, no cell is added or removed, coverage and cell count are unchanged.
Only *which* of the user's own cells carry the mass changes. No transport field, no donor, no
novelty mass, no drift shape of any kind.

λ is bisected to one global value so that pool top-1 matches S2-NN's, because the confound being
excluded is precisely "the swapped world is harder". Eight bisection steps gave **λ = 0.9941**,
top-1 **0.4966** against S2-NN's 0.4943 — matched to 0.0023, inside the pre-registered 0.01
tolerance, so the comparison is at equal difficulty.

## 2. The primary result

| condition | n | top-1 | transported-calibrator ECE |
|---|---|---|---|
| REAL (no perturbation) | 24,517 | 0.7128 | 0.0123 |
| **S4-ACC** (structureless, accuracy-matched) | 24,517 | **0.4966** | **0.0383** |
| **S2-NN** (transport swap, from O3) | 6,224 | 0.4943 | **0.1609** |

**Ratio 0.238.** Pre-registered rule: ≥0.8 → withdraw the causal reading; ≤0.5 → it survives;
between → inconclusive. **It survives.**

The reading this licenses, and no more than this: the collapse of a transported confidence
calibrator under the transport swap is **not** explained by the swapped world merely being
harder. Destroying the matcher's input until accuracy falls to the same 0.49 costs only
0.0383 of calibration error, a quarter of what the drift-shaped swap costs.

**Integrity check: passed exactly.** `max|d_app − d_own| = 0.000e+00` across all 24,517 users, so
the construction does what it claims and measured drift really is untouched.

## 3. The λ sweep

Reported because the verdict should not rest on one tuned point.

| λ | top-1 | transported ECE |
|---|---|---|
| 0.2 | 0.7195 | 0.0253 |
| 0.4 | 0.7230 | 0.0445 |
| 0.6 | 0.7141 | 0.0674 |
| 0.8 | 0.6719 | **0.0767** |
| 1.0 | 0.4826 | 0.0410 |

ECE is **non-monotone in λ** and peaks at 0.8, not at 1.0. At λ = 1.0 the profile is pure noise,
every score collapses, and the Platt features push the predicted probability down far enough that
it accidentally tracks the low accuracy again. Worth stating in the paper, because it shows that
"more disruption" and "worse calibration" are not the same axis.

## 4. The two qualifications

**4.1 The margin is a factor of 2, not 4, if you take the worst case.** The pre-registered
comparison is at matched accuracy and gives 0.0383/0.1609 = 0.238. But if the question is "could
*any* structureless perturbation of this family reach S2's figure", the relevant number is the
sweep maximum, 0.0767, which is **0.477** of 0.1609 — under the 0.5 threshold, but only just
(threshold 0.0804). The separation between "drift-shaped" and "structureless" is therefore about
twofold at the worst point of the ladder, not fourfold. The paper should quote both.

**4.2 The pre-registered secondary check failed, and the check was wrong.** I wrote:

> regressing the calibration residual on `d_app` within S4 must give approximately zero,
> **because `d_app` is constant by construction**. A non-zero coefficient means the construction
> is broken and the primary result is void.

The result is **b(d_app) = −0.0413**, with monotone per-quintile gaps +0.060, +0.043, −0.005,
−0.037, −0.064. By the letter of that rule the primary result is void.

**The rule's stated reason is false, and that is my error.** `d_app` is not constant across users
in S4; it *equals each user's own real drift*, which varies widely — the quintile means run from
0.519 to 0.923. What is constant is the *identity* `d_app ≡ d_own`, not the value. So regressing
the residual on `d_app` inside S4 is regressing it on each user's own real drift, and real drift
predicting a transported calibrator's failure is the paper's main finding (O1 §4b). A coefficient
of −0.0413 is what *must* happen if the paper is right.

The thing that rule was meant to protect against — the construction accidentally injecting
drift — is tested directly and exactly, at `max|d_app − d_own| = 0.000e+00`. That check passes.

I am recording this rather than deleting the rule, because "the check failed, so I have decided
the check was wrong" is exactly the move that should be viewed with suspicion. The distinction
that makes it defensible is that the rule's **premise** is checkable and false, and the property
it intended to verify is **separately verified and holds**. A reader who disagrees can still read
§2 and the integrity line and reach their own conclusion.

## 5. What this does not establish

- **One perturbation family.** Dirichlet mixing over a user's own support is structureless, but it
  is one kind of structureless. A perturbation that preserves the mass distribution and permutes
  *which* cell is which, or one that displaces the anchor without changing drift, would be
  different controls. S4 excludes "any large perturbation suffices"; it does not exclude every
  conceivable alternative mechanism.
- **The O3 pre-registered verdict is still INCONCLUSIVE.** S4 strengthens the *amended,
  exploratory* reading in O3 §7. It does not convert the pre-registered result, which was
  triggered by S1 and S2 disagreeing in direction, into support.
- **Arm A only.** No control was run on the matched-control arm, where the transported calibrator
  does not fail in the first place, so there is nothing to control for.
- **Different split from O3.** This script fits its own calibration/evaluation split, so its REAL
  transported ECE is 0.0123 where O3 recorded 0.0086. The comparison that matters is S4 against
  S2 *within* their own runs, and both use a transported calibrator fitted on a real world.

## 6. Reproduce

```
gzcat data/raw/yjmob100k/yjmob100k-dataset2.csv.gz | awk -f scripts/windows_profiles.awk > win_d2.tsv
python3 scripts/s4_control.py win_d2.tsv
```

Console transcript: `s4_control_output.txt`. Artefacts: `s4_control.npz`, excluded from git
because it carries per-user probabilities. The `RuntimeWarning: ... in matmul` lines are the
Accelerate BLAS artefact documented in [`ATYP_BASELINE.md`](ATYP_BASELINE.md) §8; every feature
was checked finite there.
