inline void fvmf_3d_mult_mm_geof(const DG_FP *geof, const DG_FP *mm_factor, const DG_FP *u_in, 
                                 const DG_FP *v_in, const DG_FP *w_in, DG_FP *u_out, 
                                 DG_FP *v_out, DG_FP *w_out) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  for(int m = 0; m < DG_NP; m++) {
    DG_FP u_tmp = 0.0;
    DG_FP v_tmp = 0.0;
    DG_FP w_tmp = 0.0;
    for(int n = 0; n < DG_NP; n++) {
      // int ind = m * DG_NP + n;
      int ind = DG_MAT_IND(m, n, DG_NP, DG_NP);
      u_tmp += mm_factor[n] * mass_mat[ind] * u_in[n];
      v_tmp += mm_factor[n] * mass_mat[ind] * v_in[n];
      w_tmp += mm_factor[n] * mass_mat[ind] * w_in[n];
    }
    u_out[m] += u_tmp * geof[J_IND];
    v_out[m] += v_tmp * geof[J_IND];
    w_out[m] += w_tmp * geof[J_IND];
  }
}