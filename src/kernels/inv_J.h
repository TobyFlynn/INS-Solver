inline void inv_J(const double *J, const double *tmp, double *u) {
  for(int i = 0; i < 15; i++) {
    u[i] = tmp[i] / J[i];
  }
}
