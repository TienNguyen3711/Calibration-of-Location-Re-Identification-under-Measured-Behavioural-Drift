# Premise P8 and Arm C: GeoLife

Run 2026-11-02. Closes the last open premise in [../A_PREMISES.md](../A_PREMISES.md) §4.

Scripts: [../scripts/ingest_geolife.py](../scripts/ingest_geolife.py),
[../scripts/geolife_feasibility.py](../scripts/geolife_feasibility.py),
[../scripts/armc_run.py](../scripts/armc_run.py). Output:
[P8_armC_feasibility_output.txt](P8_armC_feasibility_output.txt), [P8_armC_output.txt](P8_armC_output.txt).

Reproduce (from `Codebase/code/`):
```
python3 scripts/ingest_geolife.py "data/raw/geolife/Geolife Trajectories 1.3" geolife_bins.tsv
python3 scripts/geolife_feasibility.py geolife_bins.tsv
python3 scripts/armc_run.py geolife_bins.tsv 500 56 0 40 relative
```

**Verdict: Arm C corroborates the drift–outcome association, and cannot speak to calibration.**
Six of seven configurations give a negative `r`(drift, correct) whose 95% CI excludes zero,
including both the coarser and the finer grid that P8 demands. The one configuration that is
protocol-identical to Arms A and B is also the one with the fewest users, and there the CI contains
zero. All of that is reported below rather than the best configuration being presented alone.

## 1. Representation, fixed a priori

To match YJMob100K: **500 m cells, 30-minute slots, one observation per user per slot, the first
one chronologically** — YJMob100K's own rule. Grid origin fixed at (39.0 N, 115.0 E), equirectangular
at cos(39.9°).

**Timestamps are GMT in the `.plt` files and are shifted to Beijing local time (UTC+8) before the
day and slot are computed.** Without this the daily rhythm and any night-time definition are
displaced by eight hours. This is a correctness detail that the User Guide states and that is easy
to miss.

Ingest: 182 users, **18,669 files read** (one dropped by name: `163/Trajectory/20000101231219.plt`,
a logger that had lost its clock, per `data/raw/geolife/SOURCE.md`), **24,876,975 point lines**,
1 unparseable line, leaving **84,741 occupied 30-minute bins** — a 293× reduction.

## 2. The central representation problem: GeoLife logs trips, YJMob samples presence

| per-user quantity, after the 30-minute dedup | min | p25 | median | p75 | max |
|---|---|---|---|---|---|
| occupied bins | 2 | 43 | **128.5** | 446 | 6,987 |
| active days | 1 | 7 | 21.5 | 63 | 1,251 |
| span in days | 1 | 11 | 71.5 | 168 | 1,934 |
| **bins per active day** | 1.5 | 4.5 | **6.7** | 9.2 | 23.0 |

A YJMob100K user has 48 slots available per day and about 1,176 records over 75 days. A GeoLife
user has a median of **6.7 bins per active day** out of 48. GeoLife records journeys, not presence.
The same discretisation therefore yields a far sparser profile, and this — not the user count — is
the deepest difference between Arm C and Arms A/B. It must be stated in the paper whenever Arm C is
cited.

## 3. Usable users: computed, not assumed

A_PREMISES §2 required `L` to be chosen so that enough users keep coverage in both windows, and the
count to be computed by the build script. Two schemes were evaluated because GeoLife users are not
contemporaneous.

**CALENDAR** — one start date for everyone, as in Arms A/B. Usable users with ≥20 bins in *both*
windows, start chosen to maximise the count:

| L | g = 0 | 28 | 90 | 180 | 365 |
|---|---|---|---|---|---|
| 28 | 31 | 26 | 20 | 12 | 6 |
| 56 | 33 | 29 | 21 | 12 | 6 |
| 90 | **36** | 32 | 23 | 13 | 6 |
| 180 | 36 | 33 | 25 | 14 | 7 |

**RELATIVE** — each user's windows anchored on their own first active day. Usable users with ≥20
bins in both:

| L | g = 0 | 28 | 90 | 180 | 365 |
|---|---|---|---|---|---|
| 28 | **72** | 57 | 30 | 14 | 7 |
| 56 | 67 | 49 | 36 | 19 | 9 |
| 90 | 54 | 46 | 26 | 14 | 8 |
| 180 | 30 | 27 | 18 | 13 | 9 |

Anchoring on each user's *densest* stretch instead of their first was also tried and is **worse**
(56 at L=56, g=0 against 67), because the densest window is followed by a sparser one.

**Choice, made on coverage only, with no outcome inspected:** primary = **relative anchoring,
L = 56 days, g = 0, ≥40 bins per window → 56 users**. Rationale: the calendar scheme keeps only 36
and the pool is the binding constraint; L = 56 is eight weeks, long enough for a weekly rhythm,
where L = 28 with a ≥20-bin floor gives thin profiles. The calendar variant is reported as a
sensitivity, not hidden.

## 4. Results

Pool = the usable users themselves, so chance top-1 is 1/n. The P2 interleaved-day split is
preserved: even days feed the matcher, odd days the drift.

| cell | L | g | K<sub>min</sub> | scheme | n | top-1 | chance | `r`(JS drift, correct) | 95% CI |
|---|---|---|---|---|---|---|---|---|---|
| **500 m** | **56** | **0** | **40** | **relative** | **56** | **0.732** | 0.018 | **−0.423** | **[−0.620, −0.193]** |
| 250 m | 56 | 0 | 40 | relative | 56 | 0.750 | 0.018 | −0.378 | [−0.580, −0.140] |
| 1000 m | 56 | 0 | 40 | relative | 56 | 0.607 | 0.018 | −0.480 | [−0.672, −0.249] |
| 500 m | 90 | 0 | 20 | **calendar** | 36 | 0.778 | 0.028 | **−0.129** | **[−0.484, +0.257]** |
| 500 m | 56 | 28 | 40 | relative | 39 | 0.590 | 0.026 | −0.514 | [−0.731, −0.237] |
| 500 m | 56 | 90 | 40 | relative | 27 | 0.481 | 0.037 | −0.428 | [−0.719, −0.067] |
| 500 m | 28 | 0 | 20 | relative | 69 | 0.565 | 0.015 | −0.339 | [−0.544, −0.117] |

For comparison, tie-corrected: Arm A `r` = **−0.537**, Arm B `r` = **−0.315**. Arm C's primary
−0.423 sits between them.

**(a) P8's grid-sensitivity obligation is met and the finding survives it.** At 250 m, 500 m and
1000 m the association is −0.378, −0.423 and −0.480, every CI excluding zero. Matching accuracy does
move as expected — 0.750, 0.732, 0.607 — so a coarser grid costs linkability without destroying the
drift relationship.

**(b) The calendar configuration does not reproduce the association.** `r` = −0.129, CI
[−0.484, +0.257]. It has the fewest users (36) and the worst coverage confounding of all seven
settings (`r`(drift, coverage) = −0.403 against −0.255…−0.339 elsewhere). We do not claim it as a
failure of the hypothesis, and we do not drop it: the honest reading is that the one Arm-C
configuration that is protocol-identical to Arms A/B has too little power to decide, and the
configurations that do decide differ from Arms A/B in a stated way (non-contemporaneous pool).

**(c) The gap costs accuracy here, as in Arm B and unlike Arm A.** Top-1 falls 0.732 → 0.590 → 0.481
across g = 0, 28, 90 while the association holds at −0.423, −0.514, −0.428. Arm A's flat gap
response is the confounded one (its target window is always the shock window,
[P1_P3_P6_P7.md](P1_P3_P6_P7.md) §P3).

**(d) The P2 decoupling behaves as in Arm A.** `r`(drift on held-out days, drift on matcher days) is
**0.78 to 0.93** across configurations, so drift is a stable person-level property here too and the
interleaved split costs little. This partly satisfies acceptance criterion 1's requirement to rerun
the decoupling on Arm C.

## 5. Calibration is not reported for Arm C, and this is not a gap we can close

With a pool of 56, the Platt calibrator would be fitted on about 28 users containing roughly 20
correct matches. [P1_P3_P6_P7.md](P1_P3_P6_P7.md) §P6 showed that a **1,000**-user pool is already
degenerate for calibration in the control arm (ECE 0.354, calibrator unidentified); Arm C is two
orders of magnitude smaller. Any ECE, Brier or reliability curve computed here would be noise with a
number attached.

Consequence for the paper: **Arm C corroborates C1/C4 (drift predicts linkage outcome, under a
different drift regime and on raw GPS) and says nothing about C2, C2b or C3 (the calibration
findings).** The replication claim must be scoped to exactly that. This is a limitation of the data,
not of the design: no public location dataset we examined has both a long horizon and a pool large
enough for calibration (gap G6 in `Material/_ai_research/intake/gaps.md`).

## 6. What this does not establish

- **No calibration replication**, per §5.
- **The association is not established in the protocol-identical configuration** (§4b).
- **Beijing, 2007–2012, 182 volunteers**, mostly students and researchers — not a population sample.
  Every Arm C statement carries n.
- GeoLife has **no shock**: Arm C tests gradual, endogenous drift. It is therefore a replication of
  the Arm B regime, not of Arm A, which weakens its value as corroboration of the disruption result.
- The relative-anchoring scheme means pool candidates are **not contemporaneous**, so between-user
  seasonal and year-to-year differences enter the drift measure. Defensible for an attacker, who
  need not hold contemporaneous reference profiles, but it is a difference from Arms A/B.
- The three GeoLife citations its terms require are still **not in `references.bib`**
  (`../../paper/ETHICS.md` §2).
