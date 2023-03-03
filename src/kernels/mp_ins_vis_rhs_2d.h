inline void mp_ins_vis_rhs_2d(const DG_FP *factor, const DG_FP *rho,
                              DG_FP *vRHS0, DG_FP *vRHS1) {
  for(int i = 0; i < DG_NP; i++) {
    vRHS0[i] = (*factor) * rho[i] * vRHS0[i];
    vRHS1[i] = (*factor) * rho[i] * vRHS1[i];
  }
}
