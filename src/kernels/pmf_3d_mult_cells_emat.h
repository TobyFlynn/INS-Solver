inline void pmf_3d_mult_cells_emat(const int *order, const DG_FP *eMat,
                              const DG_FP *in0, const DG_FP *in1, const DG_FP *in2, const DG_FP *in3,
                              DG_FP *out0, DG_FP *out1, DG_FP *out2, DG_FP *out3) {
  const int p = *order;
  const DG_FP *emat_mat = &eMat[(p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  for(int i = 0; i < dg_np; i++) {
    out0[i] = 0.0;
    out1[i] = 0.0;
    out2[i] = 0.0;
    out3[i] = 0.0;
    for(int j = 0; j < DG_NUM_FACES * dg_npf; j++) {
      int ind = DG_MAT_IND(i, j, dg_np, DG_NUM_FACES * dg_npf);
      out0[i] += emat_mat[ind] * in0[j];
      out1[i] += emat_mat[ind] * in1[j];
      out2[i] += emat_mat[ind] * in2[j];
      out3[i] += emat_mat[ind] * in3[j];
    }
  }
}
