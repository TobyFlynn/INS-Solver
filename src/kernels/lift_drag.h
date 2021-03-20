inline void lift_drag(const int *bedge_type, const int *bedgeNum, const double *p,
                      const double *dQ0dx, const double *dQ0dy,
                      const double *dQ1dx, const double *dQ1dy,
                      const double *nx, const double *ny, const double *sJ,
                      double *cd, double *cl) {
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

  // If cylinder edge
  if(*bedge_type == 2) {
    for(int i = 0; i < 5; i++) {
      *cd += lift_drag_vec[i] * sJ[exInd + i] * (-p[fmask[i]] * nx[exInd + i] + nu * (nx[exInd + i] * 2.0 * dQ0dx[fmask[i]] + ny[exInd + i] * (dQ1dx[fmask[i]] + dQ0dy[fmask[i]])));
      *cl += lift_drag_vec[i] * sJ[exInd + i] * (-p[fmask[i]] * ny[exInd + i] + nu * (nx[exInd + i] * (dQ1dx[fmask[i]] + dQ0dy[fmask[i]]) + ny[exInd + i] * 2.0 * dQ1dy[fmask[i]]));
    }
  }
}
