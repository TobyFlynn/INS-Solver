inline void ins_3d_pr_3(const DG_FP *dt, const DG_FP *dpdx, const DG_FP *dpdy,
                        const DG_FP *dpdz, const DG_FP *uT, const DG_FP *vT,
                        const DG_FP *wT, DG_FP *uTT, DG_FP *vTT, DG_FP *wTT) {
  for(int i = 0; i < DG_NP; i++) {
    uTT[i] = uT[i] - *dt * dpdx[i];
    vTT[i] = vT[i] - *dt * dpdy[i];
    wTT[i] = wT[i] - *dt * dpdz[i];
  }
}
