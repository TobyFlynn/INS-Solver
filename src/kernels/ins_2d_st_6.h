inline void ins_2d_st_6(const DG_FP *rho, const DG_FP *curv, DG_FP *st_x, DG_FP *st_y) {
  for(int i = 0; i < DG_NP; i++) {
    st_x[i] *= curv[i] / (weber * rho[i]);
    st_y[i] *= curv[i] / (weber * rho[i]);
  }
}
