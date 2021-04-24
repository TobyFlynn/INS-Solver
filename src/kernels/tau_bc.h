inline void tau_bc(const int *bedgeNum, const double *J, const double *sJ, double *tau) {
  double h = 2.0 * J[FMASK[*bedgeNum * 5]] / sJ[*bedgeNum * 5];
  tau[*bedgeNum] += 100.0 * 15.0 / h;
}
