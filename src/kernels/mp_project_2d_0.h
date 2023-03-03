inline void mp_project_2d_0(const DG_FP *factor, const DG_FP *J,
                            const DG_FP *rho, const DG_FP *dpdx, const DG_FP *dpdy,
                            const DG_FP *qt0, const DG_FP *qt1,
                            DG_FP *qtt0, DG_FP *qtt1, DG_FP *rhs0,
                            DG_FP *rhs1, DG_FP *dpdn) {
  for(int i = 0; i < DG_NP; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];
    rhs0[i] = qtt0[i];
    rhs1[i] = qtt1[i];
  }

  for(int i = 0; i < DG_G_NP; i++) {
    dpdn[i] = 0.0;
  }
}
