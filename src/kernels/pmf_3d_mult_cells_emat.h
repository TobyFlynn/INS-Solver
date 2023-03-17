inline void pmf_3d_mult_cells_emat(const int *order, const DG_FP *eMat,
                              const DG_FP *in0, const DG_FP *in1, const DG_FP *in2, const DG_FP *in3,
                              DG_FP *out0, DG_FP *out1, DG_FP *out2, DG_FP *out3) {
  const int p = *order;
  const DG_FP *emat_mat = &eMat[(p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, in0, 1.0, out0);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, in1, 1.0, out1);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, in2, 1.0, out2);
  op2_in_kernel_gemv(false, dg_np, DG_NUM_FACES * dg_npf, 1.0, emat_mat, dg_np, in3, 0.0, out3);
}
