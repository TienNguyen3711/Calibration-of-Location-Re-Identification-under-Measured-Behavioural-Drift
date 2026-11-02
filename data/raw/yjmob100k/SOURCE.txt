# YJMob100K — source and verification

## Provenance
- Zenodo record `10.5281/zenodo.10836269` (v3, 18 March 2024).
- Licence: CC BY 4.0 (machine-readable `license.id = cc-by-4.0` from the Zenodo API).
- The record lists exactly the four files held here.
- Cite: `yabe2024yjmob100k` (Scientific Data 11(1), 397, 2024).

## Files and schema

| File | Schema | Rows | Users |
|---|---|---|---|
| `yjmob100k-dataset1.csv.gz` | `uid,d,t,x,y` | 111,535,175 | 100,000 |
| `yjmob100k-dataset2.csv.gz` | `uid,d,t,x,y` | 29,389,749 | 25,000 |
| `cell_POIcat.csv.gz` | `x,y,POIcategory,POI_count` | 221,159 | — |
| `POI_datacategories.csv` | one category name per line, no header | 85 | — |

- `d` = day index 0..74. `t` = 30-minute slot 0..47. `x`,`y` = 500 m grid cell, 1..200.
- Dataset 1 is 75 business-as-usual days.
- Dataset 2 is 60 business-as-usual days followed by **15 emergency days (d = 60..74)**.
- The city is not named; there are no coordinates and no real timestamps.

## Verified on 2026-11-01 (`awk` pass over every row)
- Row counts and user counts above are from the files themselves and **match the Data Descriptor
  exactly** (111,535,175 / 29,389,749 records; 100,000 / 25,000 users).
- No malformed rows (every row has 5 fields).
- No masked values: `x` or `y` equal to 999 does not occur in either file. These are complete data,
  not the challenge's masked versions.
- Value ranges are as documented: `d` 0..74, `t` 0..47, `x` 1..200, `y` 1..200 in both files.
- `POI_datacategories.csv` begins with a UTF-8 BOM; strip it when reading the first category.

## Point that constrains the experiment design

`uid` values run 0..99,999 in Dataset 1 and **0..24,999 in Dataset 2**. The ranges overlap, and
Dataset 2's ids are contiguous from zero, so the two files appear to be independently numbered:
**`uid` 7 in Dataset 1 is not the same person as `uid` 7 in Dataset 2.**

Consequence: the two datasets cannot be linked to each other. Every reference window and target
window must be taken **within** one dataset. This is not a problem for Paper 4 — Dataset 2 supplies
the shift (reference before day 60, target in days 60-74) and Dataset 1 supplies a no-shift control
over windows of the same length — but it rules out any design that treats Dataset 1 as a candidate
pool for Dataset 2 users.

The Zenodo *description* says Dataset 1 has 80,000 and Dataset 2 has 20,000 individuals, and
mentions a "90 day period". Both statements disagree with the files and with the Data Descriptor;
the files and the paper agree with each other, so the description is the one that is wrong.

## Re-download
<https://zenodo.org/records/10836269>
