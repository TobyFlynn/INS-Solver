inline void measure_zalesak_error(const DG_FP *time, const DG_FP *box,
                  const DG_FP *cub_w, const DG_FP *geof, const DG_FP *x, 
                  const DG_FP *y, const DG_FP *surface, DG_FP *l1_err, 
                  DG_FP *l2_err, DG_FP *l_max_err) {
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    if(x[i] < box[0] || x[i] > box[2] || y[i] < box[1] || y[i] > box[3]) {
      l1_err[i] = 0.0;
      l2_err[i] = 0.0;
      continue;
    }

    const DG_FP a_ls = ps2d_get_analytical_solution(*time, x[i], y[i]);
    l1_err[i] = geof[J_IND] * cub_w[i] * fabs(surface[i] - a_ls);
    l2_err[i] = geof[J_IND] * cub_w[i] * (surface[i] - a_ls) * (surface[i] - a_ls);
    *l_max_err = fmax(*l_max_err, fabs(surface[i] - a_ls));
  }
}