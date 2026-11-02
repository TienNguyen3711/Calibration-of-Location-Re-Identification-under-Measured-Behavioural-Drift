# P4 / P10 single-pass aggregates over a YJMob100K file.
# Requires rows sorted by uid,d,t (verified separately).
# Windows: REF = d 45..59, TGT = d 60..74, BAU = d 0..59.
# Home proxy = most frequent cell in night slots (t>=40 or t<16), per the Data Descriptor.
BEGIN { FS=","; REFA=45; REFB=59; TGTA=60; TGTB=74 }

function flush_user(   k,best,bc,bestn,bcn,i) {
  if (pu=="") return
  # top cell + home cell in each window
  best=""; bc=-1; for (k in rc) if (rc[k]>bc || (rc[k]==bc && k<best)) {bc=rc[k]; best=k}
  bestn=""; bcn=-1; for (k in rn) if (rn[k]>bcn || (rn[k]==bcn && k<bestn)) {bcn=rn[k]; bestn=k}
  rtop=best; rhome=bestn
  best=""; bc=-1; for (k in tc) if (tc[k]>bc || (tc[k]==bc && k<best)) {bc=tc[k]; best=k}
  bestn=""; bcn=-1; for (k in tn) if (tn[k]>bcn || (tn[k]==bcn && k<bestn)) {bcn=tn[k]; bestn=k}
  ttop=best; thome=bestn
  nref=0; for (k in rc) nref++
  ntgt=0; for (k in tc) ntgt++
  if (rrows>0) { u_ref++; sum_nref+=nref; sum_rrows+=rrows }
  if (trows>0) { u_tgt++; sum_ntgt+=ntgt; sum_trows+=trows }
  if (rrows>0 && trows>0) {
    u_both++
    if (rtop!=ttop) ch_top++
    if (rhome!="" && thome!="") { u_home_both++; if (rhome!=thome) ch_home++ }
    # Jaccard of cell sets between the two windows
    inter=0; for (k in rc) if (k in tc) inter++
    uni=nref+ntgt-inter
    if (uni>0) sum_jac+=inter/uni
  }
  split("",rc); split("",rn); split("",tc); split("",tn); rrows=0; trows=0
}

function flush_userday() {
  if (pud=="") return
  n=0; for (k in dcell) n++
  day_uc[pd_last]+=n; day_users[pd_last]++
  split("",dcell)
}

NR==1 { next }
{
  uid=$1; d=$2+0; t=$3+0; x=$4+0; y=$5+0; cell=x"_"y
  if (uid!=pu) { flush_userday(); flush_user(); pu=uid; pud="" }
  ud=uid"_"d
  if (ud!=pud) { flush_userday(); pud=ud; pd_last=d }
  dcell[cell]=1

  rows++; day_rows[d]++; slot_rows[t]++
  if (d<=59) cells_bau[cell]++
  if (d>=REFA && d<=REFB) { cells_ref[cell]++; rc[cell]++; rrows++; if (t>=40||t<16) rn[cell]++ }
  if (d>=TGTA && d<=TGTB) { cells_tgt[cell]++; tc[cell]++; trows++; if (t>=40||t<16) tn[cell]++ }
}
END {
  flush_userday(); flush_user()
  printf "ROWS\t%d\n", rows
  for (d=0; d<=74; d++) printf "DAY\t%d\t%d\t%d\t%.4f\n", d, day_rows[d]+0, day_users[d]+0, (day_users[d]>0? day_uc[d]/day_users[d] : 0)
  for (t=0; t<48; t++) printf "SLOT\t%d\t%d\n", t, slot_rows[t]+0
  printf "USERWIN\tu_ref=%d\tu_tgt=%d\tu_both=%d\tmean_cells_ref=%.3f\tmean_cells_tgt=%.3f\tmean_rows_ref=%.1f\tmean_rows_tgt=%.1f\n", \
      u_ref+0,u_tgt+0,u_both+0,(u_ref>0?sum_nref/u_ref:0),(u_tgt>0?sum_ntgt/u_tgt:0),(u_ref>0?sum_rrows/u_ref:0),(u_tgt>0?sum_trows/u_tgt:0)
  printf "CHANGE\tboth=%d\ttop_changed=%d\thome_both=%d\thome_changed=%d\tmean_jaccard=%.4f\n", \
      u_both+0, ch_top+0, u_home_both+0, ch_home+0, (u_both>0?sum_jac/u_both:0)
  for (k in cells_bau) printf "CELLBAU\t%s\t%d\n", k, cells_bau[k]
  for (k in cells_ref) printf "CELLREF\t%s\t%d\n", k, cells_ref[k]
  for (k in cells_tgt) printf "CELLTGT\t%s\t%d\n", k, cells_tgt[k]
}
