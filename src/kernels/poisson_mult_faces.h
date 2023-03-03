inline void poisson_mult_faces(const int *pL, const DG_FP *uL, const DG_FP *opL, DG_FP *rhsL,
                               const int *pR, const DG_FP *uR, const DG_FP *opR, DG_FP *rhsR) {
  const int dg_npL = DG_CONSTANTS[(*pL - 1) * DG_NUM_CONSTANTS];
  const int dg_npR = DG_CONSTANTS[(*pR - 1) * DG_NUM_CONSTANTS];

  op2_in_kernel_gemv(false, dg_npL, dg_npR, 1.0, opL, dg_npL, uR, 1.0, rhsL);
  op2_in_kernel_gemv(false, dg_npR, dg_npL, 1.0, opR, dg_npR, uL, 1.0, rhsR);
}
