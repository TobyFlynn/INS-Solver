inline void measure_l2_error_vortex(const DG_FP *time, const DG_FP *x,
                                    const DG_FP *y, const DG_FP *vel_x,
                                    DG_FP *vel_x_res) {
  const DG_FP centre_x = 0.0;
  const DG_FP centre_y = fmod(*time, (DG_FP)40.0) > 20.0 ? fmod(*time, (DG_FP)40.0) - 40.0 : fmod(*time, (DG_FP)40.0);
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP R = 1.5;
  const DG_FP S = 13.5;
  const DG_FP M = 0.4;

  for(int i = 0; i < DG_NP; i++) {
    DG_FP rel_x = x[i] - centre_x;
    DG_FP rel_y = y[i] - centre_y;
    DG_FP rel_y_1 = y[i] - (centre_y - 40.0);
    DG_FP rel_y_2 = y[i] - (centre_y + 40.0);
    if(fabs(rel_y_1) < fabs(rel_y))
      rel_y = rel_y_1;
    if(fabs(rel_y_2) < fabs(rel_y))
      rel_y = rel_y_2;

    DG_FP f = (1.0 - rel_x * rel_x - rel_y * rel_y) / (2.0 * R * R);
    DG_FP an_vel_x = (S * rel_y * exp(f)) / (2.0 * PI * R);

    if(fabs(rel_x) <= 2.0 && fabs(rel_y) <= 2.0)
      vel_x_res[i] = (vel_x[i] - an_vel_x) * (vel_x[i] - an_vel_x);
    else
      vel_x_res[i] = 0.0;
  }
}
