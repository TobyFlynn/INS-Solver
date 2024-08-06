inline void vmf_3d_mult_mm_geof(const DG_FP *geof, const DG_FP *mm_factor, const DG_FP *u,
                                const DG_FP *v, const DG_FP *w, DG_FP *u_out, DG_FP *v_out, 
                                DG_FP *w_out) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  const DG_FP J = geof[J_IND];
  for(int m = 0; m < DG_NP; m++) {
    DG_FP u_tmp = 0.0;
    DG_FP v_tmp = 0.0;
    DG_FP w_tmp = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(m, n, DG_NP, DG_NP);
      u_tmp += mass_mat[ind] * u[n];
      v_tmp += mass_mat[ind] * v[n];
      w_tmp += mass_mat[ind] * w[n];
    }
    u_out[m] += *mm_factor * u_tmp * J;
    v_out[m] += *mm_factor * v_tmp * J;
    w_out[m] += *mm_factor * w_tmp * J;
  }
}