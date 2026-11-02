"""P8: compute Arm C's usable user count instead of assuming it.

A_PREMISES §2 says L (window length in days) must be 'chosen so that enough users retain
coverage in both windows', and §4 P8 requires the discretisation to be fixed a priori with a
sensitivity check. This script supplies the numbers both decisions need.

Two window schemes are evaluated, because GeoLife users are not contemporaneous:
  CALENDAR   fixed calendar windows, as in Arms A and B. Comparable, but users must overlap
             in time, which most GeoLife users do not.
  RELATIVE   each user's windows anchored on their own first active day. Keeps far more
             users, at the cost that pool candidates are no longer contemporaneous, so
             between-user seasonal differences enter the drift. A stated difference from
             Arms A/B, not a defect -- an attacker's reference profiles need not be
             contemporaneous.

Usage: geolife_feasibility.py <geolife_bins.tsv>
"""
import sys, math, numpy as np
from collections import defaultdict

BINS = sys.argv[1]
U = defaultdict(list)
for line in open(BINS):
    uid, day, slot, lat, lon = line.split("\t")
    U[int(uid)].append((int(day), int(slot), float(lat), float(lon)))
print(f"users {len(U)}; occupied 30-min bins {sum(len(v) for v in U.values()):,}\n")

nb = np.array([len(v) for v in U.values()])
days = {u: sorted({d for d, _, _, _ in v}) for u, v in U.items()}
nd = np.array([len(days[u]) for u in U])
span = np.array([days[u][-1] - days[u][0] + 1 for u in U])
bpd = np.array([len(U[u])/len(days[u]) for u in U])

print("Per-user coverage after the 30-minute dedup")
print(f"  {'quantity':<28} {'min':>6} {'p25':>7} {'median':>7} {'p75':>7} {'max':>8}")
for nm, v in (("occupied bins", nb), ("active days", nd), ("span in days", span),
              ("bins per ACTIVE day", bpd)):
    print(f"  {nm:<28} {v.min():>6.1f} {np.quantile(v,.25):>7.1f} {np.median(v):>7.1f} "
          f"{np.quantile(v,.75):>7.1f} {v.max():>8.1f}")
print(f"\n  For comparison, YJMob100K Dataset 2 gives every user 48 possible slots a day and")
print(f"  about 1,176 records over 75 days. GeoLife's median user has {np.median(bpd):.1f} bins per")
print( "  active day: it logs TRIPS, not presence. This is the central Arm C representation issue.")

# ------------------------------------------------------------------ calendar overlap
print("\nCALENDAR scheme: how many users are simultaneously active?")
act = defaultdict(set)
for u in U:
    for d in days[u]: act[d].add(u)
allday = sorted(act)
print(f"  {'window L':>9} {'best start (day idx)':>21} {'users with >=1 bin':>19} "
      f"{'>=20 bins':>10} {'>=40 bins':>10}")
for L in (28, 56, 90, 180):
    best = None
    for s in range(allday[0], allday[-1] - 2*L, 7):
        us = set()
        for d in range(s, s + L): us |= act.get(d, set())
        if best is None or len(us) > best[1]: best = (s, len(us))
    s = best[0]
    cnt = defaultdict(int)
    for u in U:
        cnt[u] = sum(1 for d, _, _, _ in U[u] if s <= d < s + L)
    print(f"  {L:>9} {s:>21} {sum(1 for u in cnt if cnt[u] >= 1):>19} "
          f"{sum(1 for u in cnt if cnt[u] >= 20):>10} {sum(1 for u in cnt if cnt[u] >= 40):>10}")

# ------------------------------------------------------------------ relative windows
def counts(v):
    bd = defaultdict(int)
    for d, _, _, _ in v: bd[d] += 1
    return bd

BD = {u: counts(U[u]) for u in U}
def win(u, s, L):
    return sum(BD[u].get(d, 0) for d in range(s, s + L))

def usable_first(L, g, Kmin):
    k = 0
    for u in U:
        a = days[u][0]
        if win(u, a, L) >= Kmin and win(u, a + L + g, L) >= Kmin: k += 1
    return k

def usable_dense(L, g, Kmin):
    k = 0
    for u in U:
        best, bs = -1, days[u][0]
        for s in days[u]:
            c = win(u, s, L)
            if c > best: best, bs = c, s
        if best >= Kmin and win(u, bs + L + g, L) >= Kmin: k += 1
    return k

GAPS = (0, 28, 90, 180, 365)
for name, fn in (("anchored on each user's FIRST active day", usable_first),
                 ("anchored on each user's DENSEST L-day stretch", usable_dense)):
    print(f"\nRELATIVE scheme, {name}")
    print("  ref = [a, a+L), tgt = [a+L+g, a+2L+g); usable = BOTH windows hold >= Kmin bins")
    for Kmin in (20, 40):
        print(f"  Kmin = {Kmin} bins per window")
        print("    L \\ g " + "".join(f"{g:>8}" for g in GAPS))
        for L in (28, 56, 90, 180):
            print(f"    {L:>5}  " + "".join(f"{fn(L, g, Kmin):>8}" for g in GAPS))

# ------------------------------------------------------------------ calendar PAIRS
print("\n" + "="*96)
print("CALENDAR scheme, BOTH windows: ref = [s, s+L), tgt = [s+L+g, s+2L+g), same s for everyone")
print("="*96)
print("  the start s is chosen to MAXIMISE usable users, which is a design choice made on")
print("  coverage only -- no outcome is involved, so it does not select on the result")
best_cfg = {}
for Kmin in (20, 40):
    print(f"\n  Kmin = {Kmin} bins per window")
    print("    L \\\\ g " + "".join(f"{g:>11}" for g in GAPS))
    for L in (28, 56, 90, 180):
        cells = []
        for g in GAPS:
            bk, bs = -1, None
            for s in range(allday[0], allday[-1] - (2*L + g), 7):
                k = sum(1 for u in U if win(u, s, L) >= Kmin
                        and win(u, s + L + g, L) >= Kmin)
                if k > bk: bk, bs = k, s
            cells.append((bk, bs))
            best_cfg[(Kmin, L, g)] = (bk, bs)
        print(f"    {L:>5}  " + "".join(f"{k:>5}@d{s:<5}" for k, s in cells))
