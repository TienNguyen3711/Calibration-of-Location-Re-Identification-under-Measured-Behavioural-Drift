# DO NOT USE IN PAPER 4

These files are quarantined on **2026-11-01**. They must not be read by any script in this
repository, and no number derived from them may appear in Paper 4 — not as a main result, not as a
control, not as a robustness check, and not as pre-training data.

## Why

The data are the HuMob Challenge 2024 four-city release (LYMob-4Cities family). The Zenodo record
for that dataset (`10.5281/zenodo.14219563`) attaches these conditions to the CC BY 4.0 licence:

> "Researchers may use this dataset for publications and reports, as long as: 1) Users shall not
> carry out activities that involve unethical usage of the data, **including attempts at
> re-identifying data subjects**, harming individuals, or damaging companies […]"

> "Downloading this dataset implies agreement with the above two conditions."

> "This data can only be used for the purpose of participating in the Humob Challenge 2024."

Paper 4 is a re-identification study, so the first condition rules it out directly, and the purpose
limitation rules it out again.

## Why they are kept rather than deleted

They were downloaded for a different purpose and may be used for non-re-identification work later.
Deleting them is the user's decision, not an automatic one.

## Files

`city_A_alldata.csv.gz`, `city_B_alldata.csv.gz`, `city_C_alldata.csv.gz`, `city_D_alldata.csv.gz`,
`POIdata_cityA.csv`, `POIdata_cityB.csv`, `POIdata_cityC.csv`, `POIdata_cityD.csv`.

Schema `uid,d,t,x,y`, day index 1..75, timeslot 0..47, grid 1..200.
Users: 150,000 / 30,000 / 25,000 / 20,000 (A/B/C/D).

## Unresolved

The local file names and user counts do not match the Zenodo record (which states
100,000 / 25,000 / 20,000 / 6,000 and `cityA-dataset.zip` style names). The exact release these came
from is therefore unconfirmed. This does not change the conclusion, but it must be settled before
these files are used for anything.

See `../../LICENCES.md` section 3.
