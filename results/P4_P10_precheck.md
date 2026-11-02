# Pre-check of premises P4 and P10

Run 2026-11-01, before any experiment code. Tests two premises from
[../A_PREMISES.md](../A_PREMISES.md) that, had they failed, would have forced the three-arm design to
be rewritten.

**Verdict: both premises hold, and more strongly than required.** Two consequences for the design
are recorded in §4; they are not failures, but they change how drift must be measured.

## 1. Method

One streaming `awk` pass per file (`../scripts/precheck_p4_p10.awk`), relying on the rows being
sorted by `uid,d,t` (verified on the first 4M rows of each file: zero inversions). Per-user state is
flushed at each `uid` change, so memory is constant rather than proportional to the user count.
Every row of both files was read: 111,535,175 and 29,389,749.

Windows as fixed in A_PREMISES §2: `REF = d 45..59`, `TGT = d 60..74`, `BAU = d 0..59`.
Home proxy = most frequent cell in night slots (`t >= 40` or `t < 16`, i.e. 20:00-08:00), following
the Data Descriptor's own definition.

Reproduce:
```
gzcat data/raw/yjmob100k/yjmob100k-dataset1.csv.gz | awk -f scripts/precheck_p4_p10.awk > results/precheck_agg_dataset1.tsv
gzcat data/raw/yjmob100k/yjmob100k-dataset2.csv.gz | awk -f scripts/precheck_p4_p10.awk > results/precheck_agg_dataset2.tsv
python3 scripts/precheck_p4_p10_report.py results/precheck_agg_dataset1.tsv results/precheck_agg_dataset2.tsv
```
The two aggregate `.tsv` files are kept in this folder, so the report can be regenerated without
re-reading 140M rows.

## 2. P4 — Dataset 1 is a valid control for Dataset 2

> P4: Arm B is a valid control only if Dataset 1 and Dataset 2 share city, grid, resolution and
> broadly similar population behaviour, differing in the shock alone.

**Result: PASS.** The two files cover the same city **and the same calendar period**.

| Evidence | Value | Reading |
|---|---|---|
| Spatial footprint, days 0..59 in both (occupied-cell sets) | Jaccard **0.872** | same built-up area |
| Spatial footprint, days 0..59 (cell-visit proportions) | Pearson **0.972**, Spearman **0.975**, TVD 0.107 | same spatial distribution |
| Daily rhythm over the 48 half-hour slots | Pearson **0.9990**, TVD 0.0107 | same daily cycle |
| Weekly minimum of mobility | phase **0** in both (`day mod 7`) | same weekday phase |
| Correlation of the two daily mobility series, days 0..59 | **+0.990** | same calendar period |

The weekly-phase result is the decisive one. Low-mobility days fall on exactly the same day indices
in both files (0, 1, 7, 8, 14, …), so the two releases are not merely the same city — they are the
same weeks. Arm B therefore controls for city, grid, season, weekday structure and window position
simultaneously.

**Caveat to carry forward.** The Dataset 2 population is systematically **more active** than
Dataset 1: mean cells per user-day is 10-13% higher at every weekly phase, and mean records per user
in the reference window are 242.4 against 223.8. The offset is close to a constant multiplicative
factor rather than a structural difference, but it is real and must enter the models as a covariate
(premise P7), not be ignored.

## 3. P10 — the day 60..74 shock is real, and absent from Dataset 1

> P10: Arm A assumes d 60..74 of Dataset 2 is a genuine behavioural disruption. […] if the shock is
> invisible in aggregate, Arm A has no contrast with Arm B.

**Result: PASS, decisively.**

| Measure, REF (d45-59) → TGT (d60-74) | Dataset 1 (control) | Dataset 2 (shock) |
|---|---|---|
| records in window | **+3.64%** | **−9.81%** |
| mean active users per day | +2.62% | **+4.45%** |
| mean distinct cells per user-day | **+0.97%** | **−12.43%** |
| users whose top cell changed | 24,178 / 99,731 = **24.2%** | 12,473 / 24,611 = **50.7%** |
| users whose night-time home cell changed | 29,966 / 99,481 = **30.1%** | 12,264 / 24,544 = **50.0%** |
| mean Jaccard of visited-cell sets between windows | **0.2575** | **0.1981** |

The single most informative statistic is the breakdown of the between-file correlation:

```
corr(D1 daily series, D2 daily series), days  0..59 :  +0.990
corr(D1 daily series, D2 daily series), days 60..74 :  −0.336
```

The two populations move in near-perfect lockstep for sixty days and then decouple completely. The
weekly rhythm itself breaks: Dataset 2 keeps its Sunday minimum at day 56, then loses it at days 63
and 70, and instead collapses on days 66, 67 and 73 (6.90, 7.49 and 6.71 cells per user-day against
a business-as-usual level above 10). Dataset 1 keeps its regular weekly dips at 56, 63 and 70
throughout.

**The drop is not an observability artefact.** During the shock window Dataset 2 has *more* users
observed per day (+4.45%) while producing *fewer* records (−9.81%) and visiting *fewer* distinct
cells per user-day (−12.43%). Fewer people moving is therefore not the explanation; the same or more
people moved less. This partly discharges premise P7 for Arm A, though per-user coverage must still
be computed and used as a covariate.

## 4. Consequences for the design (new, not previously recorded)

**C1 — Anchor-change drift has a high noise floor.** In Dataset 1, with no shock, **30.1%** of users
change their night-time home cell between two adjacent 15-day windows, and 24.2% change their top
cell. At 500 m cells with 30-minute slots, "the person moved house" and "the estimate flipped to an
adjacent cell" are not distinguishable at the single-user level. The shock raises these to 50.0% and
50.7%, so the signal is real (about +20 points), but the anchor-change component of `δ` must be
reported with this baseline beside it, and should not be the primary drift measure.

**C2 — Baseline cell-set churn is large.** Mean Jaccard between the reference and target visited-cell
sets is only **0.258** even in the unshocked control. Of roughly 69 reference cells, about 28 are
revisited. A matcher built on exact cell-set overlap will therefore be weak on this data regardless
of drift, which argues for distribution- or histogram-based profiles (the
`tournier2022expanding` family already chosen in A_PREMISES §3) rather than set intersection.

**C3 — Arm B is a stronger control than assumed.** A_PREMISES §1 claimed B controls for "window
position". It in fact controls for city, grid, calendar period, weekday phase and window position at
once. The paper may state the contrast more strongly than originally planned.

## 5. What this does not establish

- Nothing about **per-person** drift predicting re-identification. These are population aggregates.
  The claim in A_PREMISES §0 is untouched by this check.
- Nothing about premise **P2** (drift and matcher score must not be mechanically coupled), which
  remains the highest-risk open item.
- The nature of the emergency is still undisclosed; this check shows only that *a* disruption
  occurred, not what it was.
- Dataset 1 vs Dataset 2 **individual identity** is still unknown and unusable: both files number
  users from zero (`data/raw/yjmob100k/SOURCE.md`). Nothing here changes that.
