inline void pmf_2d_mult_faces(const int *order, const int *faceNum,
                              const bool *reverse, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in,
                              const DG_FP **in_x, const DG_FP **in_y,
                              DG_FP **l_x, DG_FP **l_y, DG_FP **out) {
  const int dg_np  = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS + 1];

  const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(*order - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  const DG_FP gtau = 2.0 * (*order + 1) * (*order + 2) * fmax(fscale[0], fscale[1]);

  const bool rev = *reverse;
  for(int j = 0; j < dg_npf; j++) {
    const int fmaskL_ind = fmaskL[j];
    const int fmaskR_ind = rev ? fmaskR[dg_npf - j - 1] : fmaskR[j];
    const DG_FP diffL_u = in[0][fmaskL_ind] - in[1][fmaskR_ind];
    const DG_FP diff_u_x = in_x[1][fmaskR_ind] + in_x[0][fmaskL_ind];
    const DG_FP diff_u_y = in_y[1][fmaskR_ind] + in_y[0][fmaskL_ind];
    const DG_FP diffL_u_grad = nx[0] * diff_u_x + ny[0] * diff_u_y;
    const DG_FP diffR_u_grad = nx[1] * diff_u_x + ny[1] * diff_u_y;

    const int indL = faceNum[0] * dg_npf + j;
    out[0][indL] = 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] = nx[0] * l_tmpL;
    l_y[0][indL] = ny[0] * l_tmpL;

    const int indR = rev ? faceNum[1] * dg_npf + dg_npf - j - 1 : faceNum[1] * dg_npf + j;
    out[1][indR] = 0.5 * sJ[1] * (gtau * -diffL_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * diffL_u;
    l_x[1][indR] = nx[1] * l_tmpR;
    l_y[1][indR] = ny[1] * l_tmpR;
  }
}
