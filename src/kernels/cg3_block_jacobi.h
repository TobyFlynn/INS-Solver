inline void cg3_block_jacobi(const DG_FP *u_pre_mat, const DG_FP *v_pre_mat, const DG_FP *w_pre_mat,
                             const DG_FP *u_in, const DG_FP *v_in, const DG_FP *w_in,
                             DG_FP *u_out, DG_FP *v_out, DG_FP *w_out) {
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, u_pre_mat, DG_NP, u_in, 0.0, u_out);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, v_pre_mat, DG_NP, v_in, 0.0, v_out);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, w_pre_mat, DG_NP, w_in, 0.0, w_out);
}