inline void ins_3d_st_5(const DG_FP *curv, const DG_FP *rho, DG_FP *stx, DG_FP *sty, DG_FP *stz) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    stx[i] = stx[i] * curv[i] / (weber * rho[i]);
    sty[i] = sty[i] * curv[i] / (weber * rho[i]);
    stz[i] = stz[i] * curv[i] / (weber * rho[i]);
  }
}
