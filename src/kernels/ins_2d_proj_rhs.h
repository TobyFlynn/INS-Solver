inline void ins_2d_proj_rhs(const DG_FP *dt, const DG_FP *dpdx, const DG_FP *dpdy,
                            const DG_FP *uT, const DG_FP *vT, DG_FP *uTT,
                            DG_FP *vTT, DG_FP *proj0, DG_FP *proj1) {
  for(int i = 0; i < DG_NP; i++) {
    uTT[i] = uT[i] - *dt * dpdx[i];
    vTT[i] = vT[i] - *dt * dpdy[i];
    proj0[i] = uTT[i];
    proj1[i] = vTT[i];
  }
}
