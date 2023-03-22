inline void petsc_pre_inv_mass_dat(const DG_FP *matrix, const DG_FP *J,
                                   const DG_FP *factor, const DG_FP *in,
                                   DG_FP *out) {
  // Get constants
  // const int dg_np   = DG_CONSTANTS[(*p - 1) * 5];
  const DG_FP *mat = &matrix[(DG_ORDER - 1) * DG_NP * DG_NP];

  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0 / *J, mat, DG_NP, in, 0.0, out);

  for(int i = 0; i < DG_NP; i++) {
    out[i] /= factor[i];
  }
}
