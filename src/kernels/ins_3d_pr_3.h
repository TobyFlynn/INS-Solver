inline void ins_3d_pr_3(const double *dt, const double *dpdx, const double *dpdy,
                        const double *dpdz, const double *uT, const double *vT,
                        const double *wT, double *uTT, double *vTT, double *wTT,
                        double *dPdN) {
  for(int i = 0; i < DG_NP; i++) {
    uTT[i] = uT[i] - *dt * dpdx[i];
    vTT[i] = vT[i] - *dt * dpdy[i];
    wTT[i] = wT[i] - *dt * dpdz[i];
  }

  for(int i = 0; i < 4 * DG_NPF; i++) {
    dPdN[i] = 0.0;
  }
}
