inline void calc_sol(const double *J, double *sol) {
  for(int i = 0; i < 15; i++) {
    sol[i] = sol[i] / J[i];
  }
}
