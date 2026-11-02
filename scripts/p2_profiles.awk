# Extract per-user, per-window, per-parity-half profiles from a YJMob100K file.
# Windows: R = d 45..59 (reference), T = d 60..74 (target).
# Halves:  M = even day index (matcher half), D = odd day index (drift half).
# Disjoint by construction: no record contributes to both halves.
# Requires rows sorted by uid (verified separately).
# Output, one line per (uid, win+half) that has data:
#   U <uid> <RM|RD|TM|TD> <nrec> <ndays> | <cell>:<n> ... | <slot>:<n> ...
BEGIN { FS=","; OFS="\t"; GS[1]="RM"; GS[2]="RD"; GS[3]="TM"; GS[4]="TD"; SEP="\036" }

function dump(   gi,g,k,a,line,nd) {
  for (gi=1; gi<=4; gi++) {
    g = GS[gi]
    if (!(g in grec)) continue
    nd=0; for (k in gday) { split(k,a,SEP); if (a[1]==g) nd++ }
    line = "U\t" pu "\t" g "\t" grec[g] "\t" nd "\t|"
    for (k in gcell) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gcell[k] }
    line = line "\t|"
    for (k in gslot) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gslot[k] }
    print line
  }
  split("",grec); split("",gday); split("",gcell); split("",gslot)
}

NR==1 { next }
{
  d=$2+0
  if (d<45) next
  uid=$1
  if (uid!=pu) { if (pu!="") dump(); pu=uid }
  g = (d<=59 ? "R" : "T") (d%2==0 ? "M" : "D")
  grec[g]++
  gday[g SEP d]=1
  gcell[g SEP ($4+0 "_" $5+0)]++
  gslot[g SEP ($3+0)]++
}
END { if (pu!="") dump() }
