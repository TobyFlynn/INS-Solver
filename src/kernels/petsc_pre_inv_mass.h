inline void petsc_pre_inv_mass(const DG_FP *matrix, const DG_FP *factor, 
                               const DG_FP *J, const DG_FP *in, 
                               DG_FP *out) {
  // Get constants
  // const int dg_np   = DG_CONSTANTS[(*p - 1) * 5];
  const DG_FP *mat = &matrix[(DG_ORDER - 1) * DG_NP * DG_NP];

  DG_FP tmp[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    out[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      // int ind = i + j * DG_NP;
      int ind = DG_MAT_IND(i, j, DG_NP, DG_NP);
      out[i] += mat[ind] * in[j];
    }
    out[i] = *factor * out[i] / *J;
  }
}