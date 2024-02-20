inline void ls_driver_get_err(const DG_FP *time, const DG_FP *err_dist, const DG_FP *x, const DG_FP *y, 
                              const DG_FP *s, DG_FP *err, DG_FP *l_max) {
  const DG_FP t = *time;
  const DG_FP _err_dist = *err_dist;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP a_sol = ps2d_get_analytical_solution(t, x[i], y[i]);
    if(fabs(a_sol) <= _err_dist)
      err[i] = fabs(s[i] - a_sol);
    else
      err[i] = 0.0;

    if(err[i] > *l_max) *l_max = err[i];
  }
}