inline void fpmf_grad_3d(const int *p, const DG_FP *dr, const DG_FP *ds,
                    const DG_FP *dt, const DG_FP *u, const DG_FP *fact,
                    const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                    const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                    const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                    DG_FP *ux, DG_FP *uy, DG_FP *uz) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  DG_FP tmp_r[DG_NP], tmp_s[DG_NP], tmp_t[DG_NP];
  op2_in_kernel_gemv(false, dg_np, dg_np, 1.0, dr_mat, dg_np, u, 0.0, tmp_r);
  op2_in_kernel_gemv(false, dg_np, dg_np, 1.0, ds_mat, dg_np, u, 0.0, tmp_s);
  op2_in_kernel_gemv(false, dg_np, dg_np, 1.0, dt_mat, dg_np, u, 0.0, tmp_t);

  for(int m = 0; m < dg_np; m++) {
    ux[m] = fact[m] * (*rx * tmp_r[m] + *sx * tmp_s[m] + *tx * tmp_t[m]);
    uy[m] = fact[m] * (*ry * tmp_r[m] + *sy * tmp_s[m] + *ty * tmp_t[m]);
    uz[m] = fact[m] * (*rz * tmp_r[m] + *sz * tmp_s[m] + *tz * tmp_t[m]);
  }
}
