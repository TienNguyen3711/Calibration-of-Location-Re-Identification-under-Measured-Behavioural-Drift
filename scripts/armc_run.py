"""Arm C: run the Arm A/B protocol on GeoLife (premise P8).

Representation fixed a priori to match YJMob100K: 500 m cells, 30-minute slots, one
observation per user per slot, Beijing local time. Sensitivity at 250 m and 1000 m, as P8
requires. Grid origin is fixed at (39.0 N, 115.0 E) with an equirectangular projection at
cos(39.9 deg), so a cell is 500 m x 500 m to within the latitude range of the data.

Windows, chosen on COVERAGE ONLY (no outcome was inspected when choosing them), from
scripts/geolife_feasibility.py:
  PRIMARY   relative to each user's first active day, L = 56 days, g = 0, >= 40 bins per
            window -> 56 candidate users. Chosen because the calendar-aligned scheme keeps
            only 36 and the matcher pool is the binding constraint here.
  SENSITIVITY  calendar-aligned, L = 90, g = 0, >= 20 bins -> 36 users; and the gap sweep.
The P2 interleaved-day split is preserved: even days feed the matcher, odd days the drift.

Usage: armc_run.py <geolife_bins.tsv> [cell_m] [L] [g] [Kmin] [scheme]
"""
import sys, os, math, numpy as np
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matchlib as M

BINS = sys.argv[1]
CELL = float(sys.argv[2]) if len(sys.argv) > 2 else 500.0
L    = int(sys.argv[3]) if len(sys.argv) > 3 else 56
G    = int(sys.argv[4]) if len(sys.argv) > 4 else 0
KMIN = int(sys.argv[5]) if len(sys.argv) > 5 else 40
SCHEME = sys.argv[6] if len(sys.argv) > 6 else "relative"
LAT0, LON0, LATREF = 39.0, 115.0, 39.9
DLAT = CELL/111320.0
DLON = CELL/(111320.0*math.cos(math.radians(LATREF)))

raw = defaultdict(list)
for line in open(BINS):
    uid, day, slot, lat, lon = line.split("\t")
    raw[int(uid)].append((int(day), int(slot), float(lat), float(lon)))
for u in raw: raw[u].sort()

def cell(lat, lon):
    return f"{int((lon-LON0)/DLON)}_{int((lat-LAT0)/DLAT)}"

days = {u: sorted({d for d, _, _, _ in v}) for u, v in raw.items()}
bd = {u: defaultdict(int) for u in raw}
for u, v in raw.items():
    for d, _, _, _ in v: bd[u][d] += 1
def nwin(u, s, n): return sum(bd[u].get(d, 0) for d in range(s, s+n))

if SCHEME == "relative":
    starts = {u: days[u][0] for u in raw}
else:
    bk, bs = -1, days[min(days)][0]
    alld = sorted({d for u in raw for d in days[u]})
    for s in range(alld[0], alld[-1]-(2*L+G), 7):
        k = sum(1 for u in raw if nwin(u, s, L) >= KMIN and nwin(u, s+L+G, L) >= KMIN)
        if k > bk: bk, bs = k, s
    starts = {u: bs for u in raw}

def group(u, lo, hi, half):
    c = defaultdict(int); s = np.zeros(48); n = 0; dd = set()
    for d, sl, lat, lon in raw[u]:
        if not (lo <= d < hi): continue
        if (d % 2 == 0) != (half == "M"): continue
        c[cell(lat, lon)] += 1; s[sl] += 1; n += 1; dd.add(d)
    return dict(n=n, d=len(dd), n27=0, c=dict(c), s=s)

U = {}
for u in raw:
    a = starts[u]
    if nwin(u, a, L) < KMIN or nwin(u, a+L+G, L) < KMIN: continue
    g = {}
    for half in "MD":
        g["R"+half] = group(u, a, a+L, half)
        g["T"+half] = group(u, a+L+G, a+2*L+G, half)
    if min(x["n"] for x in g.values()) < 5: continue
    U[str(u)] = g
users = sorted(U, key=int)
N = len(users)
print(f"=== Arm C (GeoLife) | cell {CELL:.0f} m | L {L} d | g {G} d | Kmin {KMIN} | {SCHEME} ===")
print(f"usable users: {N} of 182  (both windows >= {KMIN} bins AND >= 5 bins in each of the "
      f"four window x half groups)")
if N < 12:
    print("fewer than 12 users: nothing is reportable at this setting."); sys.exit(0)
for g in ("RM", "RD", "TM", "TD"):
    v = np.array([U[u][g]["n"] for u in users])
    nc = np.array([len(U[u][g]["c"]) for u in users])
    print(f"  {g}: bins/user median {np.median(v):>6.1f} (min {v.min()}, max {v.max()});  "
          f"distinct cells median {np.median(nc):>5.1f}")

res = M.match(U, users, "RM", "TM", K=min(8, N-1))
dec = M.drift(U, users, "RD", "TD")
cou = M.drift(U, users, "RM", "TM")
cov = np.array([U[u]["RD"]["n"] + U[u]["TD"]["n"] for u in users], float)
y = res["y"]
s1, s2 = res["topk"][:, 0], res["topk"][:, 1]
cmar = (s1 - s2)/np.maximum(s1, 1e-12)
print(f"\nmatcher: top-1 {y.mean():.4f} ({int(y.sum())}/{N})  "
      f"top-3 {(res['rank']<=3).mean():.4f}  median rank {np.median(res['rank']):.0f}  "
      f"(pool = {N}, so chance top-1 = {1/N:.4f})")

rng = np.random.default_rng(20261102)
def boot_spear(a, b, reps=4000):
    v = np.array([M.spear(a[i], b[i]) for i in
                  (rng.integers(0, len(a), len(a)) for _ in range(reps))])
    v = v[np.isfinite(v)]
    return np.quantile(v, [.025, .975])
print(f"\ndrift vs outcome, with bootstrap 95% CI (n = {N}); Arm A/B values beside them")
print(f"  {'component':<11} {'r(dec,y)':>9} {'95% CI':>18} {'r(dec,c)':>9} {'r(dec,cou)':>11} "
      f"{'armA r(dec,y)':>14} {'armB':>7}")
# Arms A and B, tie-corrected (midrank) Spearman -- the 2026-11-02 first pass used an
# argsort implementation that ignored ties and attenuated every r(., y) by 0.10-0.19.
REF_A = dict(js_cells=-0.537, dissim=-0.434, anchor=-0.222, novelty=-0.436,
             d_entropy=-0.128, js_slots=-0.163)
REF_B = dict(js_cells=-0.315, dissim=None, anchor=None, novelty=None,
             d_entropy=None, js_slots=None)
for k in M.COMPS:
    lo, hi = boot_spear(dec[k], y)
    rb = REF_B.get(k)
    print(f"  {k:<11} {M.spear(dec[k],y):>+9.3f} [{lo:>+7.3f},{hi:>+7.3f}] "
          f"{M.spear(dec[k],cmar):>+9.3f} {M.spear(dec[k],cou[k]):>+11.3f} "
          f"{REF_A[k]:>+14.3f} {('%+.3f'%rb) if rb is not None else '-':>7}")
print(f"\n  r(drift, coverage) = {M.spear(dec['js_cells'],cov):+.3f}   "
      f"mean JS drift = {dec['js_cells'].mean():.3f} "
      f"(Arm A 0.702, Arm B 0.593)")
print(f"\nCALIBRATION IS NOT REPORTED for Arm C: with a pool of {N} the Platt model would be fitted")
print(f"  on about {N//2} users containing roughly {int(y.mean()*N//2)} correct matches. P6 showed a")
print( "  1,000-user pool is already degenerate for calibration; this is two orders smaller.")
