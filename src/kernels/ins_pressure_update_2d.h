inline void ins_pressure_update_2d(const DG_FP *factor, const DG_FP *dpdx,
                                   const DG_FP *dpdy, const DG_FP *qt0,
                                   const DG_FP *qt1, DG_FP *qtt0,
                                   DG_FP *qtt1, DG_FP *dpdn) {
  for(int i = 0; i < DG_NP; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i];
  }

  for(int i = 0; i < DG_G_NP; i++) {
    dpdn[i] = 0.0;
  }
}
