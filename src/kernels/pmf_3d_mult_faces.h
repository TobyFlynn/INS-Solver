inline void pmf_3d_mult_faces(const int **order, const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in, const DG_FP **in_x,
                              const DG_FP **in_y, const DG_FP **in_z, DG_FP **l_x,
                              DG_FP **l_y, DG_FP **l_z, DG_FP **out) {
  const int p = order[0][0];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const int findL = faceNum[0] * dg_npf;
  const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < dg_npf; j++) {
    const DG_FP diffL_u = in[0][fmaskL[j]] - in[1][fmaskR_corrected[j]];
    const DG_FP diffL_u_x = nx[0] * (in_x[1][fmaskR_corrected[j]] + in_x[0][fmaskL[j]]);
    const DG_FP diffL_u_y = ny[0] * (in_y[1][fmaskR_corrected[j]] + in_y[0][fmaskL[j]]);
    const DG_FP diffL_u_z = nz[0] * (in_z[1][fmaskR_corrected[j]] + in_z[0][fmaskL[j]]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[0][indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] += nx[0] * l_tmpL;
    l_y[0][indL] += ny[0] * l_tmpL;
    l_z[0][indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[1][fmaskR[j]] - in[0][fmaskL_corrected[j]];
    const DG_FP diffR_u_x = nx[1] * (in_x[1][fmaskR[j]] + in_x[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_y = ny[1] * (in_y[1][fmaskR[j]] + in_y[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_z = nz[1] * (in_z[1][fmaskR[j]] + in_z[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    const int indR = findR + j;
    out[1][indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[1][indR] += nx[1] * l_tmpR;
    l_y[1][indR] += ny[1] * l_tmpR;
    l_z[1][indR] += nz[1] * l_tmpR;
  }
}
