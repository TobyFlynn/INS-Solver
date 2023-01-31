inline void project_setup(const double *dr_, const double *ds_,
                          const double *mass_, const double *rx,
                          const double *sx, const double *ry, const double *sy,
                          double *op_0, double *op_1, double *op_2, double *op_3) {
  // Get matrices
  const double *Dr   = &dr_[(DG_ORDER - 1) * DG_NP * DG_NP];
  const double *Ds   = &ds_[(DG_ORDER - 1) * DG_NP * DG_NP];
  const double *Mass = &mass_[(DG_ORDER - 1) * DG_NP * DG_NP];

  // Div-div term
  double Dx[DG_NP * DG_NP], Dy[DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      Dx[ind] = rx[0] * Dr[ind] + sx[0] * Ds[ind];
      Dy[ind] = ry[0] * Dr[ind] + sy[0] * Ds[ind];
    }
  }

  double Dx_t[DG_NP * DG_NP], Dy_t[DG_NP * DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int op_ind = i + j * DG_NP;
      Dx_t[op_ind] = 0.0;
      Dy_t[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        int a_ind = i + k * DG_NP;
        int b_ind = j * DG_NP + k;
        Dx_t[op_ind] += Mass[a_ind] * Dx[b_ind];
        Dy_t[op_ind] += Mass[a_ind] * Dy[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int op_ind = i + j * DG_NP;
      op_0[op_ind] = 0.0;
      op_1[op_ind] = 0.0;
      op_2[op_ind] = 0.0;
      op_3[op_ind] = 0.0;
      for(int k = 0; k < DG_NP; k++) {
        int a_ind = i * DG_NP + k;
        int b_ind = j * DG_NP + k;
        op_0[op_ind] += Dx[a_ind] * Dx_t[b_ind];
        op_1[op_ind] += Dy[a_ind] * Dy_t[b_ind];
        op_2[op_ind] += Dy[a_ind] * Dx_t[b_ind];
        op_3[op_ind] += Dx[a_ind] * Dy_t[b_ind];
      }
    }
  }
}
