inline void euler_2d_wQ_2(const double *dt, double *q0, double *q1, double *q2, 
                          double *q3, const double *rk00, const double *rk01, 
                          const double *rk02, const double *rk03, const double *rk10, 
                          const double *rk11, const double *rk12, const double *rk13, 
                          const double *rk20, const double *rk21, const double *rk22, 
                          const double *rk23) {
  for(int i = 0; i < DG_NP; i++) {
    q0[i] = q0[i] + *dt * (rk00[i] / 6.0 + rk10[i] / 6.0 + rk20[i] * (2.0 / 3.0));
    q1[i] = q1[i] + *dt * (rk01[i] / 6.0 + rk11[i] / 6.0 + rk21[i] * (2.0 / 3.0));
    q2[i] = q2[i] + *dt * (rk02[i] / 6.0 + rk12[i] / 6.0 + rk22[i] * (2.0 / 3.0));
    q3[i] = q3[i] + *dt * (rk03[i] / 6.0 + rk13[i] / 6.0 + rk23[i] * (2.0 / 3.0));
  }
}