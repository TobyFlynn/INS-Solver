inline void measure_l2_error_ls_advec(const DG_FP *time, const DG_FP *x,
                                      const DG_FP *y, const DG_FP *surface,
                                      DG_FP *err) {
  const DG_FP vel_x = 1.0 / sqrt(2.0);
  const DG_FP vel_y = 1.0 / sqrt(2.0);
  const DG_FP centre_x = -0.65 + *time  * vel_x;
  const DG_FP centre_y = -0.65 + *time  * vel_y;

  for(int i = 0; i < DG_NP; i++) {
    const DG_FP a_ls = sqrt((x[i] - centre_x) * (x[i] - centre_x) + (y[i] - centre_y) * (y[i] - centre_y)) - 0.2;

    // Will always be true for now
    if(fabs(a_ls) <= 10.0)
      err[i] = (surface[i] - a_ls) * (surface[i] - a_ls);
    else
      err[i] = 0.0;
  }
}
