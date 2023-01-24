inline void euler_wQ_1(const double *dt, const double *q0, const double *q1, 
                       const double *q2, const double *q3, const double *rk00, 
                       const double *rk01, const double *rk02, const double *rk03,
                       const double *rk10, const double *rk11, const double *rk12, 
                       const double *rk13, double *wq0, double *wq1, double *wq2, 
                       double *wq3) {
  for(int i = 0; i < DG_NP; i++) {
    wq0[i] = q0[i] + *dt * (rk00[i] / 4.0 + rk10[i] / 4.0);
    wq1[i] = q1[i] + *dt * (rk01[i] / 4.0 + rk11[i] / 4.0);
    wq2[i] = q2[i] + *dt * (rk02[i] / 4.0 + rk12[i] / 4.0);
    wq3[i] = q3[i] + *dt * (rk03[i] / 4.0 + rk13[i] / 4.0);
  }
}