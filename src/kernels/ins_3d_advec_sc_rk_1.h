inline void ins_3d_advec_sc_rk_1(const DG_FP *dt, const DG_FP *uT, const DG_FP *vT,
                                 const DG_FP *wT, DG_FP *u0, DG_FP *v0, DG_FP *w0,
                                 const DG_FP *u1, const DG_FP *v1, const DG_FP *w1,
                                 DG_FP *u_out, DG_FP *v_out, DG_FP *w_out) {
  const DG_FP time_int = 0.5 * *dt / 4.0;
  for(int i = 0; i < DG_NP; i++) {
    u_out[i] = uT[i] - time_int * (u0[i] + u1[i]);
    v_out[i] = vT[i] - time_int * (v0[i] + v1[i]);
    w_out[i] = wT[i] - time_int * (w0[i] + w1[i]);
    u0[i] = (u0[i] + u1[i]) / 6.0;
    v0[i] = (v0[i] + v1[i]) / 6.0;
    w0[i] = (w0[i] + w1[i]) / 6.0;
  }
}
