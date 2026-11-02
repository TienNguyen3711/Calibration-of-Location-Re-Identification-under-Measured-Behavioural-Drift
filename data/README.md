# Paper 4 — data

Raw datasets, their provenance and their terms of use. Audited and reorganised on **2026-11-01**.

**Read `LICENCES.md` before using anything here.** One of the three datasets must not be used in
this paper at all.

## Layout

```
data/
├── README.md       this file
├── LICENCES.md     terms of use per dataset, and the decision for Paper 4
├── MANIFEST.csv    one row per file: sha256, size, rows, users, value ranges
├── raw/            untouched downloads, never written to by code
│   ├── yjmob100k/                      USABLE   506 MB
│   ├── geolife/                        USABLE   1.6 GB
│   └── _QUARANTINE_humob2024_4cities/  NOT USABLE  517 MB
├── interim/        intermediate build artefacts (empty; created by scripts)
└── processed/      analysis-ready tables (empty; created by scripts)
```

`raw/` is read-only by convention: every script reads from `raw/` and writes only to `interim/` or
`processed/`. Nothing in `raw/` is ever modified in place, so the manifest checksums stay valid.

## The three datasets

| Dataset | What it gives Paper 4 | Size | Usable |
|---|---|---|---|
| **YJMob100K** | Large n, 75 days, 500 m / 30 min grid. Dataset 2 has a **labelled shift** (15 emergency days, d = 60..74); Dataset 1 is a matched **no-shift control**. | 100,000 + 25,000 users | **yes** (CC BY 4.0) |
| **GeoLife 1.3** | Small n, multi-year, raw GPS. The **natural long-term drift** case. | 182 users, ~40-80 usable | **yes** (non-commercial) |
| **HuMob 2024 4 cities** | — | 225,000 users | **no** — terms forbid re-identification |

## What the audit established

Every row of every CSV was read (`awk`) on 2026-11-01.

**YJMob100K is clean and matches its Data Descriptor exactly.** 111,535,175 rows / 100,000 users
(Dataset 1) and 29,389,749 rows / 25,000 users (Dataset 2); no malformed rows; no masked values;
ranges `d` 0..74, `t` 0..47, `x`,`y` 1..200.

**GeoLife is complete**: 182 users, 18,670 `.plt` files, matching the User Guide. One file has a
broken date (user 163, dated 2000-01-01) and 69 users have transport-mode labels where the guide
claims 73.

Three findings change how the experiments must be built:

1. **YJMob100K Dataset 1 and Dataset 2 cannot be linked to each other.** Both number their users
   from zero (0..99,999 and 0..24,999), so the same `uid` in the two files is not the same person.
   Reference and target windows must both come from inside a single dataset.
   Details in `raw/yjmob100k/SOURCE.md`.
2. **GeoLife's median user has 19 active days.** Only 41 users span 180 days or more. The usable
   sample for a reference/gap/target design is roughly 40-80 users and must be computed, not assumed.
   Details in `raw/geolife/SOURCE.md`.
3. **The four-city files are off-limits**, and their provenance does not match the Zenodo record
   they most resemble. Details in `raw/_QUARANTINE_humob2024_4cities/DO_NOT_USE.md`.

## Known inconsistencies in the upstream documentation

Recorded because they will otherwise cause confusion later:

- The YJMob100K **Zenodo description** says Dataset 1 has 80,000 users, Dataset 2 has 20,000, and
  mentions a "90 day period". The **files and the published paper** both say 100,000 / 25,000 and
  75 days. The description is wrong.
- The LYMob-4Cities Zenodo record states 100,000 / 25,000 / 20,000 / 6,000 users for cities A-D;
  the local files hold 150,000 / 30,000 / 25,000 / 20,000, and index days 1..75 rather than 0..74.
- The GeoLife User Guide says 73 users have labels; 69 `labels.txt` files exist.

## Git

The data themselves are **not** in git (`.gitignore`). The documentation in this folder is, so
provenance, licences and checksums are versioned even though the bytes are not.

To restore this folder on another machine: download the two sources with the links in
`raw/yjmob100k/SOURCE.md` and `raw/geolife/SOURCE.md`, place them at the paths in `MANIFEST.csv`,
and check the sha256 values.
