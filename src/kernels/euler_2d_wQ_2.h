inline void euler_2d_wQ_2(const DG_FP *dt, DG_FP *q0, DG_FP *q1, DG_FP *q2, 
                          DG_FP *q3, const DG_FP *rk00, const DG_FP *rk01, 
                          const DG_FP *rk02, const DG_FP *rk03, const DG_FP *rk10, 
                          const DG_FP *rk11, const DG_FP *rk12, const DG_FP *rk13, 
                          const DG_FP *rk20, const DG_FP *rk21, const DG_FP *rk22, 
                          const DG_FP *rk23) {
  for(int i = 0; i < DG_NP; i++) {
    q0[i] = q0[i] + *dt * (rk00[i] / 6.0 + rk10[i] / 6.0 + rk20[i] * (2.0 / 3.0));
    q1[i] = q1[i] + *dt * (rk01[i] / 6.0 + rk11[i] / 6.0 + rk21[i] * (2.0 / 3.0));
    q2[i] = q2[i] + *dt * (rk02[i] / 6.0 + rk12[i] / 6.0 + rk22[i] * (2.0 / 3.0));
    q3[i] = q3[i] + *dt * (rk03[i] / 6.0 + rk13[i] / 6.0 + rk23[i] * (2.0 / 3.0));
  }
}