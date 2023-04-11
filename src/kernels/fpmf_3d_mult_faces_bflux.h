inline void fpmf_3d_mult_faces_bflux(const int *p, const int *faceL,
                          const int *faceNums, const int *fmaskL_corrected,
                          const int *fmaskR_corrected, const DG_FP *nx,
                          const DG_FP *ny, const DG_FP *nz, const DG_FP *tau,
                          const DG_FP *sJ, const DG_FP **in, const DG_FP *factor,
                          const DG_FP **in_x, const DG_FP **in_y,
                          const DG_FP **in_z, DG_FP *l_x, DG_FP *l_y, DG_FP *l_z,
                          DG_FP *out) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];
  const int faceInd = *faceL ? 0 : 1;

  const int findL = faceNums[faceInd] * dg_npf;
  const int *fmaskL = &fmask[faceNums[faceInd] * dg_npf];
  const int *fmaskR = *faceL ? fmaskR_corrected : fmaskL_corrected;

  const DG_FP int_fact = 0.5 * sJ[faceInd];

  for(int j = 0; j < dg_npf; j++) {
    const int fmaskIndL = fmaskL[j];
    const int fmaskIndR = fmaskR[j];
    const DG_FP diffL_u = in[0][fmaskIndL] - in[1][fmaskIndR];
    const DG_FP diffL_u_x = nx[faceInd] * (in_x[1][fmaskIndR] + in_x[0][fmaskIndL]);
    const DG_FP diffL_u_y = ny[faceInd] * (in_y[1][fmaskIndR] + in_y[0][fmaskIndL]);
    const DG_FP diffL_u_z = nz[faceInd] * (in_z[1][fmaskIndR] + in_z[0][fmaskIndL]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[indL] += int_fact * (*tau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = int_fact * factor[fmaskIndL] * -diffL_u;
    l_x[indL] += nx[faceInd] * l_tmpL;
    l_y[indL] += ny[faceInd] * l_tmpL;
    l_z[indL] += nz[faceInd] * l_tmpL;
  }
}
