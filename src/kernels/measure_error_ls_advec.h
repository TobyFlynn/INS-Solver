inline void measure_error_ls_advec(const DG_FP *time, const DG_FP *err_width,
                                   const DG_FP *cub_w, const DG_FP *x, const DG_FP *y, 
                                   const DG_FP *surface, DG_FP *l1_err, 
                                   DG_FP *l2_err, DG_FP *l_max_err) {
  const DG_FP vel_x = 1.0 / sqrt(2.0);
  const DG_FP vel_y = 1.0 / sqrt(2.0);
  const DG_FP centre_x = -0.65 + *time  * vel_x;
  const DG_FP centre_y = -0.65 + *time  * vel_y;
  const DG_FP _err_width = *err_width;

  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    const DG_FP a_ls = sqrt((x[i] - centre_x) * (x[i] - centre_x) + (y[i] - centre_y) * (y[i] - centre_y)) - 0.2;

    if(fabs(a_ls) <= _err_width) {
      l1_err[i] = cub_w[i] * fabs(surface[i] - a_ls);
      l2_err[i] = cub_w[i] * (surface[i] - a_ls) * (surface[i] - a_ls);
      *l_max_err = fmax(*l_max_err, fabs(surface[i] - a_ls));
    } else {
      l1_err[i] = 0.0;
      l2_err[i] = 0.0;
    }
  }
}
