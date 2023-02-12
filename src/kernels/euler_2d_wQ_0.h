inline void euler_2d_wQ_0(const double *dt, const double *q0, const double *q1, 
                          const double *q2, const double *q3, const double *rk00, 
                          const double *rk01, const double *rk02, const double *rk03,
                          double *wq0, double *wq1, double *wq2, double *wq3) {
  for(int i = 0; i < DG_NP; i++) {
    wq0[i] = q0[i] + *dt * rk00[i];
    wq1[i] = q1[i] + *dt * rk01[i];
    wq2[i] = q2[i] + *dt * rk02[i];
    wq3[i] = q3[i] + *dt * rk03[i];
  }
}