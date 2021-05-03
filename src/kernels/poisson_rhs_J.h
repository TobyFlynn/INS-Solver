inline void poisson_rhs_J(const double *J, double *qx, double *qy,
                          double *gradx, double *grady) {
  for(int i = 0; i < 15; i++) {
    qx[i] = qx[i] / J[i];
    qy[i] = qy[i] / J[i];
    gradx[i] = gradx[i] / J[i];
    grady[i] = grady[i] / J[i];
  }
}
