# Calibration of Location Re-Identification under Measured Behavioural Drift

Code and written record for a privacy-measurement study: a location re-identification
matcher is given a probabilistic confidence output, and the reliability of that confidence
is measured **conditional on how much each individual's behaviour changed** between a
reference and a target period.

The study is a measurement of how much a risk estimate can be trusted. It is not an attempt
to re-identify any person, and it produces no per-person output.

## What is here

```
A_PREMISES.md        the design contract: 11 premises, each tested or bounded, with a change log
B_SWAP_DESIGN.md     pre-registration of the drift-swap counterfactual, with its declared amendment
scripts/             one script per analysis, plus matchlib.py, the shared matcher/drift/calibration library
results/*.md         the written record, one file per premise or objective
results/*_output.txt the console transcript each record was written from
data/                provenance, licences and checksums only -- no data
```

## Design in one paragraph

Three arms. **Arm A** is the primary: a city-scale release whose final 15 days fall in a
declared emergency period, so a real behavioural disruption is present. **Arm B** is a
matched control from the same city, grid, calendar period and weekday phase, with no
disruption. **Arm C** is a replication on a second, independent source. Profiles are built
over five 15-day windows. The matcher is a deliberately standard cell-visit histogram with a
Hellinger transform, IDF weighting and cosine similarity -- it is instrumentation, not a
contribution. Its confidence is produced two ways: an exact full-pool softmax over a
temperature grid, and a Platt logistic fit on score order statistics. Drift is measured from
**odd-indexed days while the matcher consumes only even-indexed days**, so no record feeds
both; this split is what makes the drift-confidence relationship something other than
circular, and its cost is quantified in `results/P2_decoupling.md`.

## Reproducing

Each `results/*.md` ends with the exact commands that regenerate its artefacts. The raw data
is **not** redistributed; `data/raw/*/SOURCE.md` records where each source comes from and
under which terms. The profile extraction is a single streaming `awk` pass per file over tens
of millions of rows and takes a few minutes.

Python 3.9+, numpy and scipy. No network access is required for any analysis.

## Data use and exclusions

- One source is **quarantined**: its terms forbid re-identification attempts, so it is
  unused and no script reads its path. See `data/raw/_QUARANTINE_humob2024_4cities/DO_NOT_USE.md`.
- **No raw data is in this repository.**
- **No per-user artefacts are in this repository.** `results/*.npz` and `results/*.pkl` hold
  per-user re-identification outcomes, confidences, drift values and rarity scores, and are
  excluded by `.gitignore` on ethics grounds, not merely for size. `results/precheck_agg_*.tsv`
  are aggregates but contain absolute cell coordinates that together form the spatial density
  signature of a city the release does not name, and are excluded for the same reason.
  `.gitignore` states the full rationale.
- The ethics position, including open author actions and residual risks that cannot be
  discharged, is in `../paper/ETHICS.md`.

## Status

All eleven premises in `A_PREMISES.md` are resolved. The drift-swap counterfactual returned
`INCONCLUSIVE` on its pre-registered reading; the amendment that followed is labelled
exploratory and its missing control is named as missing in `results/O3_swap.md`. Three
experiments are still owed and are listed in `../paper/OUTLINE.md`: an equal-magnitude
perturbation control for the swap, drift-aware recalibration, and a second matcher family.

Negative and inconclusive results are recorded here rather than dropped, and so are
corrections to earlier runs -- including a tie-correction to every rank correlation, which
changed several reported coefficients and left every calibration metric unchanged.
