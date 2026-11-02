# Licences and terms of use

Checked on **2026-11-01** against the authoritative sources named below (Zenodo REST API for the
Zenodo records, so the text is the record's own, not a search summary). Re-check before submission.

This file records what was found and **what it means for Paper 4**, whose subject is
*re-identification of individuals in location data*. That purpose makes the "re-identification"
clauses below decisive, so each dataset is marked USABLE or NOT USABLE for this paper.

---

## 1. YJMob100K — **USABLE**

- Record: Zenodo `10.5281/zenodo.10836269` (v3, 18 March 2024), title
  "YJMob100K: City-Scale and Longitudinal Dataset of Anonymized Human Mobility Trajectories".
- Machine-readable licence field returned by the Zenodo API: `{"id": "cc-by-4.0"}`
  → **Creative Commons Attribution 4.0 International**.
- Files in the record: `yjmob100k-dataset1.csv.gz`, `yjmob100k-dataset2.csv.gz`,
  `cell_POIcat.csv.gz`, `POI_datacategories.csv` — exactly the four files held in
  `raw/yjmob100k/`.
- The record description contains **no** ethical-use clause, **no** prohibition on
  re-identification, and **no** purpose limitation. CC BY 4.0 permits use and redistribution
  with attribution.
- Required citation (Data Descriptor):
  Yabe, T., Tsubouchi, K., Shimizu, T., Sekimoto, Y., Sezaki, K., Moro, E., & Pentland, A. (2024).
  *YJMob100K: City-scale and longitudinal dataset of anonymized human mobility trajectories.*
  Scientific Data, 11(1), 397. <https://doi.org/10.1038/s41597-024-03237-9>
  (bib key in this project: `yabe2024yjmob100k`)

### Caveat to address in the paper (not a licence restriction)

The Data Descriptor states, in its *Limitations of the data* section:

> "The YJMob100K dataset is intended to serve the scientific community as a benchmark data for
> human mobility prediction tasks, not for urban science or human behavioral studies."

This is the authors' statement of **intended use**, not a licence term; the licence is CC BY 4.0,
which does not restrict the field of use. The context of the sentence is that the anonymisation
(no city name, no latitude/longitude, no real POI names) makes the data unfit for urban science or
behavioural inference, not that such work is forbidden.

Paper 4 is neither a prediction study nor an urban-science study: it is a privacy evaluation of the
released data. Two points should still be stated explicitly in the paper's ethics section, because a
reviewer will raise them:

1. Paper 4 measures per-person behavioural **change** as a covariate. That is instrumentally a
   behavioural quantity, so the paper should say why this is a privacy measurement rather than a
   behavioural study of the population.
2. The re-identification is **within the released dataset only** (pseudonymous `uid` to pseudonymous
   `uid`). The data carry no city, no coordinates and no timestamps, so no result can be carried back
   to a natural person. This should be stated as the reason the work does not create new risk for
   LY Corporation's users.

A courtesy e-mail to the corresponding authors (T. Yabe, K. Tsubouchi) before submission is
advisable, and cheap. Not yet done.

---

## 2. GeoLife Trajectories 1.3 — **USABLE** (terms partly unresolved)

- Source: Microsoft Research Asia, GeoLife project; Microsoft download page id 52367.
- **No licence file is present in the downloaded archive**, and `User Guide-1.3.pdf` (3 pages,
  version 1.3, 2012/08/01) contains no licence section — only Data Description, What's new,
  Paper Citation, Data Format and Contact.
- The Microsoft download page itself displays no licence text.
- Secondary indication only: *User Guide 1.2* reportedly carried a "Microsoft Research License
  Agreement" headed "Non-Commercial Use Only". **This was not verified first-hand** and the 1.2
  guide is not in this repository.
- The User Guide 1.3 names the intended research fields, and privacy is among them:

  > "This trajectory dataset can be used in many research fields, such as mobility pattern mining,
  > user activity recognition, location-based social networks, **location privacy**, and location
  > recommendation."

  This is a direct statement by the data providers that **location-privacy research is an intended
  use**, which covers Paper 4.
- Required citation — the guide asks for three papers:
  1. Zheng, Zhang, Xie, Ma. *Mining interesting locations and travel sequences from GPS trajectories.* WWW 2009, 791-800.
  2. Zheng, Li, Chen, Xie, Ma. *Understanding mobility based on GPS data.* UbiComp 2008, 312-321.
  3. Zheng, Xie, Ma. *GeoLife: A collaborative social networking service among user, location and trajectory.* IEEE Data Engineering Bulletin 33(2), 2010, 32-40.
  None of these three is yet in `Codebase/paper/references.bib`; they must be added before writing.

### Practical position
Academic, non-commercial use with the three citations is consistent with everything found.
**Do not redistribute the raw data** (it is already excluded by `.gitignore`), since the
redistribution terms could not be confirmed.

---

## 3. HuMob Challenge 2024 / LYMob-4Cities (`city_A`-`city_D`) — **NOT USABLE for Paper 4**

Held in `raw/_QUARANTINE_humob2024_4cities/`. Not deleted, because the files may be legitimately
used for a different, non-re-identification purpose later.

- Nearest matching record: Zenodo `10.5281/zenodo.14219563`, "LYMob-4Cities: Multi-City Human
  Mobility Dataset", licence field `cc-by-4.0`, files `cityA-dataset.zip` … `POIdata_cityD.zip`.
- The record description adds conditions on top of the licence. Verbatim:

  > "Researchers may use this dataset for publications and reports, as long as: 1) **Users shall not
  > carry out activities that involve unethical usage of the data, including attempts at
  > re-identifying data subjects**, harming individuals, or damaging companies, and 2) The Data
  > Descriptor paper of an earlier version of the dataset (citation below) needs to be cited when
  > using the data for research and/or commercial purposes."

  > "Downloading this dataset implies agreement with the above two conditions."

  > "This data can only be used for the purpose of participating in the Humob Challenge 2024."

- Paper 4 is a re-identification study. The first condition prohibits exactly that, and the last
  sentence limits the purpose to the 2024 challenge. **These files must not be used in Paper 4**,
  including as a control, a pre-training corpus, or a robustness check.

### Open point on provenance
The local files do not match that record exactly, so the precise origin is unconfirmed:

| | Zenodo LYMob-4Cities says | Local files contain |
|---|---|---|
| file names | `cityA-dataset.zip` … | `city_A_alldata.csv.gz` … |
| users per city (A/B/C/D) | 100,000 / 25,000 / 20,000 / 6,000 | 150,000 / 30,000 / 25,000 / 20,000 |
| day index | (not stated) | 1..75 (YJMob100K uses 0..74) |

The local copies are probably the challenge's own "all data" release (unmasked answer keys:
no `x`/`y` value of 999 appears anywhere). Whichever release they came from, the HuMob Challenge 2024
terms quoted above are the ones the Zenodo record attaches to this data, so the conclusion is
unchanged. If these files are ever needed, their exact origin and terms must be established first.

---

## Summary

| Dataset | Location | Licence | Re-identification allowed | Paper 4 |
|---|---|---|---|---|
| YJMob100K | `raw/yjmob100k/` | CC BY 4.0 | not restricted | **USABLE** |
| GeoLife 1.3 | `raw/geolife/` | MSR, non-commercial (unconfirmed) | privacy named as intended use | **USABLE** |
| HuMob 2024 4-cities | `raw/_QUARANTINE_.../` | CC BY 4.0 **+ extra conditions** | **explicitly prohibited** | **NOT USABLE** |

## To do before submission
- [ ] Add the three GeoLife citations to `Codebase/paper/references.bib`.
- [ ] Optional courtesy e-mail to the YJMob100K corresponding authors.
- [ ] Re-check both Zenodo records (records can be versioned and re-licensed).
- [ ] Write the ethics section covering the two points in section 1 above.
