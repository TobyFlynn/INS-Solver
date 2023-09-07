inline void ins_2d_proj_5(const DG_FP *_pen, const DG_FP *geof, const DG_FP *u,
                          const DG_FP *v, DG_FP *out0, DG_FP *out1) {
  const DG_FP pen = *_pen;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP r = out0[i]; DG_FP s = out1[i];
    out0[i] = geof[RX_IND] * r + geof[SX_IND] * s;
    out1[i] = geof[RY_IND] * r + geof[SY_IND] * s;
    out0[i] *= pen;
    out1[i] *= pen;
  }

  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  op2_in_kernel_gemv(false, DG_NP, DG_NP, geof[J_IND], mass_mat, DG_NP, u, 1.0, out0);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, geof[J_IND], mass_mat, DG_NP, v, 1.0, out1);
}
