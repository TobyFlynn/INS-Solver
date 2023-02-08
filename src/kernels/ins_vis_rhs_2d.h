inline void ins_vis_rhs_2d(const double *factor, double *vRHS0, double *vRHS1) {
  for(int i = 0; i < DG_NP; i++) {
    vRHS0[i] = (*factor) * vRHS0[i];
    vRHS1[i] = (*factor) * vRHS1[i];
  }
}
