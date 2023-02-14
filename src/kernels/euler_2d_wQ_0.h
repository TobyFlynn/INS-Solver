inline void euler_2d_wQ_0(const DG_FP *dt, const DG_FP *q0, const DG_FP *q1, 
                          const DG_FP *q2, const DG_FP *q3, const DG_FP *rk00, 
                          const DG_FP *rk01, const DG_FP *rk02, const DG_FP *rk03,
                          DG_FP *wq0, DG_FP *wq1, DG_FP *wq2, DG_FP *wq3) {
  for(int i = 0; i < DG_NP; i++) {
    wq0[i] = q0[i] + *dt * rk00[i];
    wq1[i] = q1[i] + *dt * rk01[i];
    wq2[i] = q2[i] + *dt * rk02[i];
    wq3[i] = q3[i] + *dt * rk03[i];
  }
}