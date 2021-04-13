inline void poisson_rhs_qbc(const int *bedge_type, const int *bedgeNum,
                           const int *neumann0, const int *neumann1, const int *neumann2,
                           const double *u, double *fluxU) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedge_type == *neumann0 || *bedge_type == *neumann1 || *bedge_type == *neumann2) {
    // Do nothing, numerical flux should be 0 for dirichlet BCs
  } else {
    // So du is 0 if no dirichlet BCs on this edge
    for(int i = 0; i < 7; i++) {
      fluxU[exInd + i] += u[exInd + i];
    }
  }
}
