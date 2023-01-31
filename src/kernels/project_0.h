inline void project_0(const double *factor, const double *J,
                      const double *rho, const double *dpdx,
                      const double *dpdy, const double *qt0,
                      const double *qt1, double *qtt0, double *qtt1,
                      double *rhs0, double *rhs1, double *dpdn, double *prBC) {
  for(int i = 0; i < DG_NP; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];
    rhs0[i] = qtt0[i];
    rhs1[i] = qtt1[i];
  }

  for(int i = 0; i < DG_G_NP; i++) {
    prBC[i] = 0.0;
    dpdn[i] = 0.0;
  }
}
