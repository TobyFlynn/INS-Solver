inline void diff_mult(const double *vis, double *dx, double *dy) {
  double k = sqrt(*vis);
  for(int i = 0; i < DG_NP; i++) {
    dx[i] *= k;
    dy[i] *= k;
  }
}