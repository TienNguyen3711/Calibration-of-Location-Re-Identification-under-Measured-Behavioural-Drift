"""P4 / P10 analysis from the two aggregate files."""
import sys, numpy as np

def load(path):
    day, slot = {}, {}
    cells = {"BAU": {}, "REF": {}, "TGT": {}}
    userwin, change, rows = {}, {}, None
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if p[0] == "ROWS": rows = int(p[1])
        elif p[0] == "DAY": day[int(p[1])] = (int(p[2]), int(p[3]), float(p[4]))
        elif p[0] == "SLOT": slot[int(p[1])] = int(p[2])
        elif p[0].startswith("CELL"): cells[p[0][4:]][p[1]] = int(p[2])
        elif p[0] == "USERWIN": userwin = dict(kv.split("=") for kv in p[1:])
        elif p[0] == "CHANGE": change = dict(kv.split("=") for kv in p[1:])
    return dict(rows=rows, day=day, slot=slot, cells=cells, userwin=userwin, change=change)

def vec(c, keys):
    v = np.array([c.get(k, 0) for k in keys], float)
    return v / v.sum()

def compare_maps(a, b, label):
    keys = sorted(set(a) | set(b))
    pa, pb = vec(a, keys), vec(b, keys)
    occ_a, occ_b = set(a), set(b)
    jac = len(occ_a & occ_b) / len(occ_a | occ_b)
    pear = np.corrcoef(pa, pb)[0, 1]
    # Spearman without scipy
    ra = np.argsort(np.argsort(pa)); rb = np.argsort(np.argsort(pb))
    spear = np.corrcoef(ra, rb)[0, 1]
    tvd = 0.5 * np.abs(pa - pb).sum()
    cos = float(pa @ pb / (np.linalg.norm(pa) * np.linalg.norm(pb)))
    print(f"  {label:<34} cells_a={len(occ_a):6d} cells_b={len(occ_b):6d} "
          f"Jaccard={jac:.3f} Pearson={pear:.3f} Spearman={spear:.3f} cosine={cos:.3f} TVD={tvd:.3f}")
    return dict(jac=jac, pear=pear, spear=spear, tvd=tvd, cos=cos)

d1 = load(sys.argv[1]); d2 = load(sys.argv[2])
print(f"rows: dataset1={d1['rows']:,}  dataset2={d2['rows']:,}\n")

print("=" * 100)
print("P4 - are Dataset 1 and Dataset 2 the same city / comparable population?")
print("=" * 100)
print(" spatial footprint, business-as-usual days 0..59 only (same day range in both):")
compare_maps(d1["cells"]["BAU"], d2["cells"]["BAU"], "D1 BAU vs D2 BAU")
print("\n within-dataset sanity (REF vs TGT inside the SAME dataset):")
compare_maps(d1["cells"]["REF"], d1["cells"]["TGT"], "D1 ref(45-59) vs D1 tgt(60-74)")
compare_maps(d2["cells"]["REF"], d2["cells"]["TGT"], "D2 ref(45-59) vs D2 tgt(60-74)")

s1 = np.array([d1["slot"][t] for t in range(48)], float); s1 /= s1.sum()
s2 = np.array([d2["slot"][t] for t in range(48)], float); s2 /= s2.sum()
print(f"\n daily rhythm over 48 slots: Pearson={np.corrcoef(s1,s2)[0,1]:.4f}  "
      f"TVD={0.5*np.abs(s1-s2).sum():.4f}  (1.0 / 0.0 = identical)")

print("\n per-user activity in the two windows:")
for name, d in (("Dataset 1", d1), ("Dataset 2", d2)):
    u = d["userwin"]
    print(f"  {name}: users_ref={u['u_ref']:>6} users_tgt={u['u_tgt']:>6} both={u['u_both']:>6} "
          f"| mean cells ref={float(u['mean_cells_ref']):6.2f} tgt={float(u['mean_cells_tgt']):6.2f} "
          f"| mean rows ref={float(u['mean_rows_ref']):7.1f} tgt={float(u['mean_rows_tgt']):7.1f}")

print()
print("=" * 100)
print("P10 - is the day 60..74 shock visible in Dataset 2, and absent in Dataset 1?")
print("=" * 100)

def win_stats(d, a, b):
    rows = sum(d["day"][x][0] for x in range(a, b + 1))
    users = np.mean([d["day"][x][1] for x in range(a, b + 1)])
    uc = np.mean([d["day"][x][2] for x in range(a, b + 1)])
    return rows, users, uc

for name, d in (("Dataset 1 (no shock expected)", d1), ("Dataset 2 (shock at d60-74)", d2)):
    r_ref, u_ref, c_ref = win_stats(d, 45, 59)
    r_tgt, u_tgt, c_tgt = win_stats(d, 60, 74)
    print(f"\n {name}")
    print(f"   rows/window         ref={r_ref:>12,}  tgt={r_tgt:>12,}   change={100*(r_tgt-r_ref)/r_ref:+6.2f}%")
    print(f"   mean active users/d ref={u_ref:>12,.0f}  tgt={u_tgt:>12,.0f}   change={100*(u_tgt-u_ref)/u_ref:+6.2f}%")
    print(f"   mean cells/user-day ref={c_ref:>12.3f}  tgt={c_tgt:>12.3f}   change={100*(c_tgt-c_ref)/c_ref:+6.2f}%")
    ch = d["change"]; uw = d["userwin"]
    nb = int(ch["both"]); hb = int(ch["home_both"])
    print(f"   users in both windows                 : {nb:,}")
    print(f"   top cell changed ref->tgt             : {int(ch['top_changed']):,} ({100*int(ch['top_changed'])/nb:.1f}%)")
    print(f"   home cell (night) changed ref->tgt    : {int(ch['home_changed']):,} ({100*int(ch['home_changed'])/hb:.1f}% of {hb:,})")
    print(f"   mean Jaccard of visited-cell sets     : {float(ch['mean_jaccard']):.4f}")

print("\n per-day series, days 55..74 (rows | active users | mean cells per user-day):")
print("   day |        D1 rows   D1 users  D1 cells |        D2 rows   D2 users  D2 cells")
for x in range(55, 75):
    a, b = d1["day"][x], d2["day"][x]
    mark = " <<" if x == 60 else ""
    print(f"   {x:>3} | {a[0]:>13,} {a[1]:>10,} {a[2]:>9.3f} | {b[0]:>13,} {b[1]:>10,} {b[2]:>9.3f}{mark}")

print()
print("=" * 100)
print("Calendar alignment - do the two files cover the same weekly phase / period?")
print("=" * 100)
a = np.array([d1["day"][i][2] for i in range(75)])
b = np.array([d2["day"][i][2] for i in range(75)])
for name, v in (("D1", a), ("D2", b)):
    m = [v[i::7][:9].mean() for i in range(7)]
    print(f"  {name} mean cells/user-day by (day mod 7), BAU days only: "
          + " ".join(f"{x:5.2f}" for x in m) + f"   -> weekly minimum at phase {int(np.argmin(m))}")
for lbl, sl in (("days 0..59 (both business-as-usual)", slice(0, 60)), ("days 60..74 (shock window)", slice(60, 75))):
    print(f"  corr(D1 daily series, D2 daily series) over {lbl:<36}: {np.corrcoef(a[sl], b[sl])[0,1]:+.3f}")
