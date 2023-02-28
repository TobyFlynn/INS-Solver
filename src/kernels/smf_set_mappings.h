inline void smf_set_mappings(const int **order, const int **glb_ind,
                             const int *faceNum, const int *fmaskL_corrected,
                             const int *fmaskR_corrected, int *mappings) {
  const int p = order[0][0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const int findL = faceNum[0] * dg_npf;
  const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  for(int i = 0; i < dg_npf; i++) {
    mappings[4 * i + 0] = glb_ind[0][0] * DG_NP + fmaskL[i];
    mappings[4 * i + 1] = glb_ind[1][0] * DG_NP + fmaskR_corrected[i];
    mappings[4 * i + 2] = glb_ind[0][0] * DG_NUM_FACES * DG_NPF + findL + i;
    for(int j = 0; j < dg_npf; j++) {
      if(fmaskR[j] == fmaskR_corrected[i])
        mappings[4 * i + 3] = glb_ind[1][0] * DG_NUM_FACES * DG_NPF + findR + j;
    }
  }
}
