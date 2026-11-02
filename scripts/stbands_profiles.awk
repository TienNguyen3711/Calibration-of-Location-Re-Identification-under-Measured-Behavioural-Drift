# Spatio-temporal profiles: cell x time-band, for the SECOND MATCHER FAMILY.
#
# The primary matcher profiles a user by WHICH cells they visit and discards WHEN. Two people
# who visit the same places at different hours are indistinguishable to it. This pass emits the
# joint distribution instead: the key is "<x>_<y>@<band>", with four six-hour bands over the
# 48 half-hour slots.
#
#   band 0  slots  0..11   00:00-06:00
#   band 1  slots 12..23   06:00-12:00
#   band 2  slots 24..35   12:00-18:00
#   band 3  slots 36..47   18:00-24:00
#
# Everything else is held fixed so the comparison attributes to the REPRESENTATION alone: the
# same five windows, the same even/odd parity split (P2), the same downstream scoring rule.
#
# Output format is identical to windows_profiles.awk, so matchlib.parse_windows reads it
# unchanged; the cell field simply carries composite keys.
#   U <uid> <W{0..4}{M|D}> <nrec> <ndays> <n27> | <cell@band>:<n> ... | <slot>:<n> ...
BEGIN { FS=","; OFS="\t"; SEP="\036" }

function dump(   g,k,a,line,nd,seen) {
  for (g in gseen) {
    nd=0; for (k in gday) { split(k,a,SEP); if (a[1]==g) nd++ }
    line = "U\t" pu "\t" g "\t" grec[g] "\t" nd "\t" (g SEP 27 in gday ? gd27[g] : 0) "\t|"
    for (k in gcell) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gcell[k] }
    line = line "\t|"
    for (k in gslot) { split(k,a,SEP); if (a[1]==g) line = line " " a[2] ":" gslot[k] }
    print line
  }
  split("",grec); split("",gday); split("",gcell); split("",gslot)
  split("",gseen); split("",gd27)
}

NR==1 { next }
{
  d=$2+0; t=$3+0
  uid=$1
  if (uid!=pu) { if (pu!="") dump(); pu=uid }
  g = "W" int(d/15) (d%2==0 ? "M" : "D")
  gseen[g]=1
  grec[g]++
  gday[g SEP d]=1
  if (d==27) gd27[g]++
  gcell[g SEP ($4+0 "_" $5+0 "@" int(t/12))]++
  gslot[g SEP t]++
}
END { if (pu!="") dump() }
