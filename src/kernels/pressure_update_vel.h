inline void pressure_update_vel(const double *factor, const double *rho,
                                const double *dpdx, const double *dpdy,
                                const double *qt0, const double *qt1,
                                double *qtt0, double *qtt1, double *dpdn,
                                double *prBC) {
  for(int i = 0; i < DG_NP; i++) {
    // qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
    // qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];
    qtt0[i] = qt0[i] - dpdx[i];
    qtt1[i] = qt1[i] - dpdy[i];
  }

  for(int i = 0; i < DG_G_NP; i++) {
    prBC[i] = 0.0;
    dpdn[i] = 0.0;
  }
}
