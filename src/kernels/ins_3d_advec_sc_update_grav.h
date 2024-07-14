inline void ins_3d_advec_sc_update_grav(const DG_FP *a0, const DG_FP *a1,
                                   const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                                   const DG_FP *u0, const DG_FP *v0, const DG_FP *w0,
                                   DG_FP *u1, DG_FP *v1, DG_FP *w1) {
  const DG_FP grav_term = *dt * (1.0 / (froude * froude)) * (*b0 + *b1);
  for(int i = 0; i < DG_NP; i++) {
    u1[i] = *a0 * u0[i] + *a1 * u1[i];
    v1[i] = *a0 * v0[i] + *a1 * v1[i] - grav_term;
    w1[i] = *a0 * w0[i] + *a1 * w1[i];
  }
}
