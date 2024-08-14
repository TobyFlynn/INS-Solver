inline void measure_error_euler_vortex(const DG_FP *time, const DG_FP *weights, const DG_FP *geof, 
                                       const DG_FP *x, const DG_FP *y, DG_FP *rho) {
  const DG_FP centre_x = 0.0;
  const DG_FP centre_y = fmod(*time, 40.0) > 20.0 ? fmod(*time, 40.0) - 40.0 : fmod(*time, 40.0);
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP R = 1.5;
  const DG_FP S = 13.5;
  const DG_FP M = 0.4;

  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    DG_FP rel_x = x[i] - centre_x;
    DG_FP rel_y = y[i] - centre_y;
    if(centre_y > 18.0 && y[i] < -18.0) {
      rel_y = y[i] - (centre_y - 40.0);
    }
    if(centre_y < -18.0 && y[i] > 18.0) {
      rel_y = y[i] - (centre_y + 40.0);
    }

    DG_FP f = (1.0 - rel_x * rel_x - rel_y * rel_y) / (2.0 * R * R);
    DG_FP rho_tmp = 1.0 - ((S * S * M * M * (1.4 - 1.0) * exp(2.0 * f)) / (8.0 * PI * PI));
    DG_FP rho_ = pow(rho_tmp, 1.0 / (1.4 - 1.0));

    if(fabs(rel_x) <= 2.0 && fabs(rel_y) <= 2.0)
      rho[i] = weights[i] * geof[J_IND] * (rho[i] - rho_) * (rho[i] - rho_);
    else
      rho[i] = 0.0;
  }
}