inline void poisson_matrix_free_3d_mult_faces_new_1(const int **order,
                              const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP **in, const DG_FP **in_x, const DG_FP **in_y,
                              const DG_FP **in_z, const DG_FP **l_x, const DG_FP **l_y,
                              const DG_FP **l_z, DG_FP **s) {
  const int p = order[0][0];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  // const int findL = faceNum[0] * dg_npf;
  // const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  for(int i = 0; i < dg_npf; i++) {
    int find = faceNum[0] * dg_npf + i;
    const DG_FP inR   = in[1][fmaskR_corrected[i]];
    const DG_FP in_xR = in_x[1][fmaskR_corrected[i]];
    const DG_FP in_yR = in_y[1][fmaskR_corrected[i]];
    const DG_FP in_zR = in_z[1][fmaskR_corrected[i]];
    const DG_FP l_xR  = l_x[1][fmaskR_corrected[i]];
    const DG_FP l_yR  = l_y[1][fmaskR_corrected[i]];
    const DG_FP l_zR  = l_z[1][fmaskR_corrected[i]];

    double flux = 0.5 * nx[0] * (in_xR - in_x[0][fmaskL[i]]);
    flux += 0.5 * ny[0] * (in_yR - in_y[0][fmaskL[i]]);
    flux += 0.5 * nz[0] * (in_zR - in_z[0][fmaskL[i]]);
    flux += gtau * (inR - in[0][fmaskL[i]]);
    // flux -= nx[0] * l_x[0][fmaskL[i]] + ny[0] * l_y[0][fmaskL[i]] + nz[0] * l_z[0][fmaskL[i]];

    s[0][find] += fscale[0] * flux;
  }

  for(int i = 0; i < dg_npf; i++) {
    int find = faceNum[1] * dg_npf + i;
    const DG_FP inL   = in[0][fmaskL_corrected[i]];
    const DG_FP in_xL = in_x[0][fmaskL_corrected[i]];
    const DG_FP in_yL = in_y[0][fmaskL_corrected[i]];
    const DG_FP in_zL = in_z[0][fmaskL_corrected[i]];
    const DG_FP l_xL  = l_x[0][fmaskL_corrected[i]];
    const DG_FP l_yL  = l_y[0][fmaskL_corrected[i]];
    const DG_FP l_zL  = l_z[0][fmaskL_corrected[i]];

    double flux = 0.5 * nx[1] * (in_xL - in_x[1][fmaskR[i]]);
    flux += 0.5 * ny[1] * (in_yL - in_y[1][fmaskR[i]]);
    flux += 0.5 * nz[1] * (in_zL - in_z[1][fmaskR[i]]);
    flux += gtau * (inL - in[1][fmaskR[i]]);
    // flux -= nx[1] * l_x[1][fmaskR[i]] + ny[1] * l_y[1][fmaskR[i]] + nz[1] * l_z[1][fmaskR[i]];

    s[1][find] += fscale[1] * flux;
  }
}
