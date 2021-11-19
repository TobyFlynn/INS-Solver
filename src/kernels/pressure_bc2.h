inline void pressure_bc2(const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, double *prBC) {
  int exInd = *bedgeNum * DG_GF_NP;

  const double PI = 3.141592653589793238463;

  if(*problem == 1) {
    if(*bedge_type == 1) {
      // Outflow
      for(int i = 0; i < DG_GF_NP; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        prBC[exInd + i] += -cos(2.0 * PI * x1) * cos(2.0 * PI * y1) * exp(-nu * 8.0 * PI * PI * *t);
      }
    }
  }
}
