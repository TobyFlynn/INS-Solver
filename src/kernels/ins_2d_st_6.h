inline void ins_2d_st_6(const DG_FP *mat, const DG_FP *curv, DG_FP *st_x, DG_FP *st_y) {
  for(int i = 0; i < DG_NP; i++) {
    // DG_FP curv_ = 1.0 / 7.5;
    // st_x[i] *= curv_ / weber;
    // st_y[i] *= curv_ / weber;
    st_x[i] *= curv[i] / weber;
    st_y[i] *= curv[i] / weber;
  }
/*
  DG_FP tmp[DG_NP];

  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, mat, DG_NP, st_x, 0.0, tmp);

  for(int i = 0; i < DG_NP; i++) {
    st_x[i] = tmp[i];
  }

  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, mat, DG_NP, st_y, 0.0, tmp);

  for(int i = 0; i < DG_NP; i++) {
    st_y[i] = tmp[i];
  }
*/
}