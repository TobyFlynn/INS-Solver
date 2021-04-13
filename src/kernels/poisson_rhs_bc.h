inline void poisson_rhs_bc(const int *bedge_type, const int *bedgeNum,
                           const int *dirichlet0, const int *dirichlet1, const int *dirichlet2,
                           const double *u, double *fluxU) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedge_type == *dirichlet0 || *bedge_type == *dirichlet1 || *bedge_type == *dirichlet2) {
    // Do nothing, numerical flux should be 0 for dirichlet BCs
  } else {
    // So du is 0 if no dirichlet BCs on this edge
    for(int i = 0; i < 7; i++) {
      fluxU[exInd + i] += u[exInd + i];
    }
  }
}
