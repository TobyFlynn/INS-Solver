inline void vmf_2d_mult_mm_geof(const int *p, const DG_FP *geof,
                                const DG_FP *mm_factor, const DG_FP *u,
                                const DG_FP *v, DG_FP *u_out, DG_FP *v_out) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  const DG_FP J = geof[J_IND];
  for(int m = 0; m < dg_np; m++) {
    DG_FP tmp = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      tmp += mass_mat[ind] * u[n];
    }
    u_out[m] += *mm_factor * tmp * J;
  }

  for(int m = 0; m < dg_np; m++) {
    DG_FP tmp = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      tmp += mass_mat[ind] * v[n];
    }
    v_out[m] += *mm_factor * tmp * J;
  }
}