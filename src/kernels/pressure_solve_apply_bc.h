inline void pressure_solve_apply_bc(const int *edgeType, const int *edgeNum,
                             const int *d0, const int *d1, const int *d2,
                             const double *mD0, const double *mD1,
                             const double *mD2, const double *sJ,
                             const double *h, const double *tau, const double *gRho, const double *rho,
                             const double *bc, double *b) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;
  // Get right matrices for this edge
  // (TODO switch matrices to be defined on edge set instead of cell set)
  const double *mD, *gVM;
  if(*edgeNum == 0) {
    mD  = mD0;
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    mD  = mD1;
    gVM = gFInterp1_g;
  } else {
    mD  = mD2;
    gVM = gFInterp2_g;
  }

  // Calculate penalty parameter
  double tauA[7];
  for(int i = 0; i < 7; i++) {
    int ind = *edgeNum  * 7 + i;
    tauA[i] = 10 * 0.5 * 5 * 6 * (*h / gRho[ind]);
  }


  double op[7 * 15];
  // First edge term
  // gVM'*gw*rho^-1*gDnM
  for(int i = 0; i < 7 * 15; i++) {
    int indT = (i % 7) * 15 + i / 7;
    int indSJ = *edgeNum * 7 + (i % 7);
    int indRho = (i / 7);

    op[i] = gVM[indT] * gaussW_g[i % 7] * sJ[indSJ] * tauA[i % 7]
            - (1.0 / rho[indRho]) * mD[indT] * gaussW_g[i % 7] * sJ[indSJ];

    // op[i] = gVM[indT] * gaussW_g[i % 7] * sJ[indSJ] * tauA[i % 7]
    //         - (1.0 / gRho[indSJ]) * mD[indT] * gaussW_g[i % 7] * sJ[indSJ];

    // op[i] = gVM[indT] * gaussW_g[i % 7] * sJ[indSJ] * tau[*edgeNum]
    //         - (1.0 / rho[indRho]) * mD[indT] * gaussW_g[i % 7] * sJ[indSJ];
  }

  // Multiply u by ops and add to rhs
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 7; j++) {
      int op_ind = i * 7 + j;
      int bc_ind = *edgeNum * 7 + j;
      b[i] += op[op_ind] * bc[bc_ind];
    }
  }
}
