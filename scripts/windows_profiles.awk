# Per-user, per-window, per-parity-half profiles over ALL FIVE 15-day windows.
# Window w = int(d/15): W0=d0..14, W1=d15..29, W2=d30..44, W3=d45..59, W4=d60..74.
# Half M = even day index (matcher), D = odd day index (drift). Disjoint by construction (P2).
# Day 27 is flagged by the Data Descriptor (p. 5); it is EMITTED but counted separately so a
# sensitivity run can drop it. Requires rows sorted by uid.
# Output: U <uid> <W#><M|D> <nrec> <ndays> <nrec_d27> | <cell>:<n> ... | <slot>:<n> ...
BEGIN { FS=","; OFS="\t"; SEP="\036"
        n=0; for (w=0; w<5; w++) { GS[++n]="W" w "M"; GS[++n]="W" w "D" } NG=n }

function dump(   gi,g,k,a,line,nd) {
  for (gi=1; gi<=NG; gi++) {
    g = GS[gi]
    if (!(g in grec)) continue
    nd=0; for (k in gday) { split(k,a,SEP); if (a[1]==g) nd++ }
    line = "U\t" pu "\t" g "\t" grec[g] "\t" nd "\t" (g in g27 ? g27[g] : 0) "\t|"
    for (k in gcell) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gcell[k] }
    line = line "\t|"
    for (k in gslot) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gslot[k] }
    print line
  }
  split("",grec); split("",gday); split("",gcell); split("",gslot); split("",g27)
}

NR==1 { next }
{
  uid=$1; d=$2+0
  if (uid!=pu) { if (pu!="") dump(); pu=uid }
  g = "W" int(d/15) (d%2==0 ? "M" : "D")
  grec[g]++
  if (d==27) g27[g]++
  gday[g SEP d]=1
  gcell[g SEP ($4+0 "_" $5+0)]++
  gslot[g SEP ($3+0)]++
}
END { if (pu!="") dump() }
