inline void poisson_rhs_bc(const int *bedge_type, const int *bedgeNum,
                           const int *dirichlet0, const int *dirichlet1,
                           const double *dBC, const double *U, double *exU) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(*bedge_type == *dirichlet0 || *bedge_type == *dirichlet1) {
    for(int i = 0; i < 5; i++) {
      // exU[exInd + i] += -U[fmask[i]] + 2.0 * dBC[exInd + i];
      // exU[exInd + i] += dBC[exInd + i];
      exU[exInd + i] += U[fmask[i]];
    }
  } else {
    // So du is 0 if no dirichlet BCs on this edge
    // TODO check this
    for(int i = 0; i < 5; i++) {
      exU[exInd + i] += -U[fmask[i]];
    }
  }
}
