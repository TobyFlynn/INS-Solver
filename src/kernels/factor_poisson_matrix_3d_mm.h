inline void factor_poisson_matrix_3d_mm(const int *order, const double *factor,
                                        const double *mass, const double *J,
                                        double *op1) {
  const double *mass_mat = &mass[(*order - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*order - 1) * 2];

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = i + j * dg_np;
      op1[op_ind] += factor[i] * J[0] * mass_mat[op_ind];
    }
  }
}
