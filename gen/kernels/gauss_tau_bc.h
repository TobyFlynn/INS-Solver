inline void gauss_tau_bc(const int *bedgeNum, const double *fscale, double *tau) {
  tau[*bedgeNum] += 20 * 25 * fscale[*bedgeNum * 4];
}
