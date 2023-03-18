inline void mp_ins_3d_advec_3(const DG_FP *a0, const DG_FP *a1, const DG_FP *b0,
                           const DG_FP *b1, const DG_FP *dt, const DG_FP *u,
                           const DG_FP *v, const DG_FP *w, const DG_FP *u_old,
                           const DG_FP *v_old, const DG_FP *w_old,
                           const DG_FP *n0, const DG_FP *n1, const DG_FP *n2,
                           const DG_FP *n0_old, const DG_FP *n1_old,
                           const DG_FP *n2_old, const DG_FP *f0, const DG_FP *f1,
                           const DG_FP *f2, const DG_FP *f0_old,
                           const DG_FP *f1_old, const DG_FP *f2_old, DG_FP *uT,
                           DG_FP *vT, DG_FP *wT) {
  for(int i = 0; i < DG_NP; i++) {
    uT[i] = (*a0 * u[i] + *a1 * u_old[i]) - *dt * (*b0 * (n0[i] + f0[i]) + *b1 * (n0_old[i] + f0_old[i]));
    vT[i] = (*a0 * v[i] + *a1 * v_old[i]) - *dt * (*b0 * (n1[i] + f1[i]) + *b1 * (n1_old[i] + f1_old[i]));
    wT[i] = (*a0 * w[i] + *a1 * w_old[i]) - *dt * (*b0 * (n2[i] + f2[i]) + *b1 * (n2_old[i] + f2_old[i]));
  }
}
