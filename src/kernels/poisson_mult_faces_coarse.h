inline void poisson_mult_faces_coarse(const DG_FP *uL, const DG_FP *opL, DG_FP *rhsL,
                                      const DG_FP *uR, const DG_FP *opR, DG_FP *rhsR) {
  op2_in_kernel_gemv(false, DG_NP_N1, DG_NP_N1, 1.0, opL, DG_NP_N1, uR, 1.0, rhsL);
  op2_in_kernel_gemv(false, DG_NP_N1, DG_NP_N1, 1.0, opR, DG_NP_N1, uL, 1.0, rhsR);
}
