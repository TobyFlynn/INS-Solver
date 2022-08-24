inline void calc_u_bar(const double *t_n, const double *t_n_1, const double *t_new,
                       const double *u_n, const double *v_n, 
                       const double *u_n_1, const double *v_n_1,
                       double *u_bar, double *v_bar) {
  const double t_diff = *t_n - *t_n_1;
  for(int i = 0; i < DG_NP; i++) {
    double grad_u = (u_n[i] - u_n_1[i]) / t_diff;
    double c_u = u_n[i] - grad_u * (*t_n); 
    u_bar[i] = grad_u * (*t_new) + c_u;

    double grad_v = (v_n[i] - v_n_1[i]) / t_diff;
    double c_v = v_n[i] - grad_v * (*t_n); 
    v_bar[i] = grad_v * (*t_new) + c_v;
  }
}