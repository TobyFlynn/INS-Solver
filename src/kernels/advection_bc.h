// TODO double check these
inline void advection_bc(const int *bedge_type, const int *bedgeNum,
                        const double *nx, const double *ny, const double *q0,
                        const double *q1, double *exQ0, double *exQ1) {
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

  if(*bedge_type == 0) {
    // Inflow
    for(int i = 0; i < 5; i++) {
      exQ0[exInd + i] += bc_u;
      exQ1[exInd + i] += bc_v;
    }
  } else if(*bedge_type == 1) {
    // Outflow
    for(int i = 0; i < 5; i++) {
      exQ0[exInd + i] += bc_u;
      exQ1[exInd + i] += bc_v;
    }
  } else {
    // Wall
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exQ0[exInd + i] += q0[qInd] - 2 * (nx[exInd + i] * q0[qInd] + ny[exInd + i] * q1[qInd]) * nx[exInd + i];
      exQ1[exInd + i] += q1[qInd] - 2 * (nx[exInd + i] * q0[qInd] + ny[exInd + i] * q1[qInd]) * ny[exInd + i];
    }
  }
}
