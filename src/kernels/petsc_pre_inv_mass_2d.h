inline void petsc_pre_inv_mass_2d(const DG_FP *factor,
                                  const DG_FP *geof, const DG_FP *in,
                                  DG_FP *out) {
  // Get constants
  // const int dg_np   = DG_CONSTANTS[(*p - 1) * 5];
  const DG_FP *mat = dg_InvMass_kernel + (DG_ORDER - 1) * DG_NP * DG_NP;

  const DG_FP mat_fact = *factor / geof[J_IND];
  op2_in_kernel_gemv(false, DG_NP, DG_NP, mat_fact, mat, DG_NP, in, 0.0, out);
}
