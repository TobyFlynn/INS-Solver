inline void ins_3d_advec_sc_rk_0(const DG_FP *dt, const DG_FP *uT, const DG_FP *vT,
                                 const DG_FP *wT, const DG_FP *u0, const DG_FP *v0,
                                 const DG_FP *w0, DG_FP *u_out, DG_FP *v_out, DG_FP *w_out) {
  for(int i = 0; i < DG_NP; i++) {
    u_out[i] = uT[i] - *dt * u0[i];
    v_out[i] = vT[i] - *dt * v0[i];
    w_out[i] = wT[i] - *dt * w0[i];
  }
}
