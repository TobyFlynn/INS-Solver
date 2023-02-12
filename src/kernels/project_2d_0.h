inline void project_2d_0(const double *factor, const double *J,
                         const double *dpdx, const double *dpdy, 
                         const double *qt0, const double *qt1, 
                         double *qtt0, double *qtt1, double *rhs0, 
                         double *rhs1, double *dpdn) {
  for(int i = 0; i < DG_NP; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i];
    rhs0[i] = qtt0[i];
    rhs1[i] = qtt1[i];
  }

  for(int i = 0; i < DG_G_NP; i++) {
    dpdn[i] = 0.0;
  }
}