inline void poisson_matrix_free_3d_mult_faces_new_0(const int **order,
                              const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP **in, DG_FP **t_x, DG_FP **t_y, DG_FP **t_z) {
  const int p = order[0][0];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  // const int findL = faceNum[0] * dg_npf;
  // const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  for(int i = 0; i < dg_npf; i++) {
    int find = faceNum[0] * dg_npf + i;

    const DG_FP flux = 0.5 * fscale[0] * (in[1][fmaskR_corrected[i]] - in[0][fmaskL[i]]);

    t_x[0][find] += nx[0] * flux;
    t_y[0][find] += ny[0] * flux;
    t_z[0][find] += nz[0] * flux;
  }

  for(int i = 0; i < dg_npf; i++) {
    int find = faceNum[1] * dg_npf + i;

    const DG_FP flux = 0.5 * fscale[1] * (in[0][fmaskL_corrected[i]] - in[1][fmaskR[i]]);

    t_x[1][find] += nx[1] * flux;
    t_y[1][find] += ny[1] * flux;
    t_z[1][find] += nz[1] * flux;
  }
}
