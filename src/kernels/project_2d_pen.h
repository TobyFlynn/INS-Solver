inline void project_2d_pen(const DG_FP *factor, const DG_FP *a0,
                           const DG_FP *a1, const DG_FP *geof, const DG_FP *u0,
                           const DG_FP *v0, const DG_FP *u1, const DG_FP *v1,
                           const DG_FP *h, DG_FP *pen) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP J = geof[J_IND];
  DG_FP tmp0[DG_NP], tmp1[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP _u = *a0 * u0[i] + *a1 * u1[i];
    const DG_FP _v = *a0 * v0[i] + *a1 * v1[i];
    tmp0[i] = sqrt(_u * _u + _v * _v);
  }
  op2_in_kernel_gemv(false, DG_NP, DG_NP, J, mass_mat, DG_NP, tmp0, 0.0, tmp1);
  DG_FP vel = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    vel += tmp1[i];
  }
  vel /= (*h) * (*h) * (*h);

  *pen = *factor * (*h) * vel;
}
