inline void ins_3d_pr_0(const double *t, const int *bc_type, const int *faceNum,
                        const double *nx, const double *ny, const double *nz,
                        const double *fscale, const double *x, const double *y,
                        const double *z, const double *n0, const double *n1,
                        const double *n2, const double *curl20,
                        const double *curl21, const double *curl22, double *dPdN) {
  // TODO handle different boundary conditions
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const double PI = 3.141592653589793238463;

  if(*bc_type == 0 || *bc_type == 2) {
    for(int i = 0; i < DG_NPF; i++) {
      double res0 = -n0[fmaskB[i]] - curl20[fmaskB[i]] / r_ynolds;
      double res1 = -n1[fmaskB[i]] - curl21[fmaskB[i]] / r_ynolds;
      double res2 = -n2[fmaskB[i]] - curl22[fmaskB[i]] / r_ynolds;

      dPdN[fInd + i] += *fscale * (*nx * res0 + *ny * res1 + *nz * res2);
    }

    if(*bc_type == 0) {
      for(int i = 0; i < DG_NPF; i++) {
        // dPdN[fInd + i] += *fscale * (PI * cos(PI * (*t)) * (y[fmask[i]] * (1.0 - y[fmask[i]]))  * (z[fmask[i]] * (1.0 - z[fmask[i]])));
        dPdN[fInd + i] += *fscale * PI * cos(PI * (*t));
      }
    }
  }
}
