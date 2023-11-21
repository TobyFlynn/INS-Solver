inline void ins_2d_st_7(const DG_FP *curv, const DG_FP *rho, DG_FP *stx, DG_FP *sty) {
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    stx[i] = stx[i] * curv[i] / (weber * rho[i]);
    sty[i] = sty[i] * curv[i] / (weber * rho[i]);
  }
}
