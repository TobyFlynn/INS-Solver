inline void ins_vis_mm_2d(const DG_FP *mm, DG_FP *rhs0, DG_FP *rhs1) {
  DG_FP tmp0[DG_NP], tmp1[DG_NP];

  for(int m = 0; m < DG_NP; m++) {
    tmp0[m] = 0.0;
    tmp1[m] = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, DG_NP, DG_NP);
      tmp0[m] += rhs0[n] * mm[ind];
      tmp1[m] += rhs1[n] * mm[ind];
    }
  }

  for(int m = 0; m < DG_NP; m++) {
    rhs0[m] = tmp0[m];
    rhs1[m] = tmp1[m];
  }
}
