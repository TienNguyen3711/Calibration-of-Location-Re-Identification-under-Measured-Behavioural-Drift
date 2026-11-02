"""Ingest GeoLife 1.3 into the same (uid, day, slot, location) form as YJMob100K.

Rules chosen to match YJMob100K so one protocol runs on both sources (premise P8):
  * one observation per user per 30-minute bin, the FIRST one chronologically
    (YJMob100K's own rule, yabe2024yjmob100k Fig. 1);
  * 48 slots per day;
  * timestamps in the .plt files are GMT, so they are shifted to Beijing local time
    (UTC+8) before the day and slot are computed -- without this the daily rhythm and any
    night-time definition are displaced by eight hours.

Output is resolution-INDEPENDENT: raw lat/lon per occupied bin, so the grid can be changed
afterwards without re-reading 18,670 files. One expensive pass, cheap re-gridding.

  uid <TAB> day <TAB> slot <TAB> lat <TAB> lon
  day = days since 2007-01-01 in Beijing local time.

Dropped, with counts reported: the one file dated outside the documented collection period
(163/Trajectory/20000101231219.plt, a logger that had lost its clock, see data/raw/geolife/
SOURCE.md) and any unparseable line.

Usage: ingest_geolife.py <geolife_root> <out.tsv>
"""
import sys, os, datetime

ROOT, OUT = sys.argv[1], sys.argv[2]
EPOCH = datetime.date(2007, 1, 1)
SHIFT = datetime.timedelta(hours=8)
BJ = (39.0, 41.2, 115.0, 117.6)          # lat_lo, lat_hi, lon_lo, lon_hi, for reporting only
DROP_FILES = {"20000101231219.plt"}

data = os.path.join(ROOT, "Data")
uids = sorted(d for d in os.listdir(data) if d.isdigit())
n_files = n_lines = n_bad = n_bins = n_out_bj = n_dropped_files = 0
fout = open(OUT, "w")
for uid in uids:
    tdir = os.path.join(data, uid, "Trajectory")
    if not os.path.isdir(tdir): continue
    seen = set()
    rows = []
    for fn in sorted(os.listdir(tdir)):
        if not fn.endswith(".plt"): continue
        if fn in DROP_FILES:
            n_dropped_files += 1
            continue
        n_files += 1
        with open(os.path.join(tdir, fn), errors="replace") as fh:
            for i, line in enumerate(fh):
                if i < 6: continue
                n_lines += 1
                p = line.rstrip("\n").split(",")
                if len(p) < 7: n_bad += 1; continue
                try:
                    lat = float(p[0]); lon = float(p[1])
                    y, mo, da = (int(x) for x in p[5].split("-"))
                    hh, mi = int(p[6][:2]), int(p[6][3:5])
                except Exception:
                    n_bad += 1; continue
                if not (-90 <= lat <= 90 and -180 <= lon <= 180): n_bad += 1; continue
                t = datetime.datetime(y, mo, da, hh, mi) + SHIFT
                day = (t.date() - EPOCH).days
                slot = (t.hour*60 + t.minute)//30
                k = (day, slot)
                if k in seen: continue
                seen.add(k)
                rows.append((day, slot, lat, lon))
                if not (BJ[0] <= lat <= BJ[1] and BJ[2] <= lon <= BJ[3]): n_out_bj += 1
    rows.sort()
    for day, slot, lat, lon in rows:
        fout.write(f"{int(uid)}\t{day}\t{slot}\t{lat:.6f}\t{lon:.6f}\n")
    n_bins += len(rows)
fout.close()
print(f"users            {len(uids)}")
print(f"files read       {n_files:,}  (dropped by name: {n_dropped_files})")
print(f"point lines      {n_lines:,}")
print(f"unparseable      {n_bad:,}")
print(f"occupied 30-min bins written  {n_bins:,}  "
      f"(compression {n_lines/max(n_bins,1):.1f}x)")
print(f"bins outside the Beijing box {BJ}: {n_out_bj:,} ({100*n_out_bj/max(n_bins,1):.2f}%)")
