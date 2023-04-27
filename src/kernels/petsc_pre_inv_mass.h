inline void petsc_pre_inv_mass(const DG_FP *factor, const DG_FP *J,
                               const DG_FP *in, DG_FP *out) {
  // Get constants
  // const int dg_np   = DG_CONSTANTS[(*p - 1) * 5];
  const DG_FP *mat = &dg_InvMass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  op2_in_kernel_gemv(false, DG_NP, DG_NP, *factor / *J, mat, DG_NP, in, 0.0, out);
}
