inline void ins_advec_sc_update_grav_2d(const DG_FP *a0, const DG_FP *a1,
                      const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                      const DG_FP *u0, const DG_FP *v0, DG_FP *u1, DG_FP *v1) {
  const DG_FP grav_term = *dt * (1.0 / (froude * froude)) * (*b0 + *b1);
  for(int i = 0; i < DG_NP; i++) {
    u1[i] = *a0 * u0[i] + *a1 * u1[i];
    v1[i] = *a0 * v0[i] + *a1 * v1[i] - grav_term;
  }
}
