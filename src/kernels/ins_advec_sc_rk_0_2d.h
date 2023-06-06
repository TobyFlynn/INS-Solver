inline void ins_advec_sc_rk_0_2d(const DG_FP *dt, const DG_FP *uT, const DG_FP *vT,
                                 const DG_FP *u0, const DG_FP *v0, DG_FP *u_out,
                                 DG_FP *v_out) {
  for(int i = 0; i < DG_NP; i++) {
    u_out[i] = uT[i] - *dt * u0[i];
    v_out[i] = vT[i] - *dt * v0[i];
  }
}
