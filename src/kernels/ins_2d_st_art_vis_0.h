inline void ins_2d_st_art_vis_0(const int *faceNum, const bool *reverse,
                               const DG_FP *nx, const DG_FP *ny,
                               const DG_FP *sJ, const DG_FP **u,
                               const DG_FP **v, DG_FP **out, int **out_count) {
  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * DG_NPF];
  const int *fmaskR = &fmask[faceNumR * DG_NPF];
  const bool rev = *reverse;

  const DG_FP *emat = &dg_Emat_kernel[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF * DG_NP];

  DG_FP maxVelU = 0.0;
  DG_FP maxVelV = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    maxVelU = fmax(maxVelU, fabs(u[0][fmaskL_ind]));
    maxVelU = fmax(maxVelU, fabs(u[1][fmaskR_ind]));
    maxVelV = fmax(maxVelV, fabs(v[0][fmaskL_ind]));
    maxVelV = fmax(maxVelV, fabs(v[1][fmaskR_ind]));
  }

  DG_FP w_L[DG_NPF], w_R[DG_NPF];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    w_L[i] = 0.0;
    for(int j = 0; j < DG_NPF; j++) {
      int f_ind = faceNumL * DG_NPF + j;
      int ind = DG_MAT_IND(fmaskL_ind, f_ind, DG_NP, DG_NPF * DG_NUM_FACES);
      w_L[i] += emat[ind];
    }

    const int fmaskR_ind = fmaskR[i];
    w_R[i] = 0.0;
    for(int j = 0; j < DG_NPF; j++) {
      int f_ind = faceNumR * DG_NPF + j;
      int ind = DG_MAT_IND(fmaskR_ind, f_ind, DG_NP, DG_NPF * DG_NUM_FACES);
      w_R[i] += emat[ind];
    }
  }

  DG_FP disSenU_L = 0.0;
  DG_FP disSenV_L = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    disSenU_L += w_L[i] * (u[0][fmaskL_ind] - u[1][fmaskR_ind]);
    disSenV_L += w_L[i] * (v[0][fmaskL_ind] - v[1][fmaskR_ind]);
  }

  DG_FP disSenU_R = 0.0;
  DG_FP disSenV_R = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    disSenU_R += w_R[i] * (u[1][fmaskR_ind] - u[0][fmaskL_ind]);
    disSenV_R += w_R[i] * (v[1][fmaskR_ind] - v[0][fmaskL_ind]);
  }

  disSenU_L = fabs(disSenU_L * sJ[0] / (maxVelU));
  disSenV_L = fabs(disSenV_L * sJ[0] / (maxVelV));
  disSenU_R = fabs(disSenU_R * sJ[1] / (maxVelU));
  disSenV_R = fabs(disSenV_R * sJ[1] / (maxVelV));

  if(maxVelU < 1e-8) {
    disSenU_L = 0.0;
    disSenU_R = 0.0;
  }

  if(maxVelV < 1e-8) {
    disSenV_L = 0.0;
    disSenV_R = 0.0;
  }

  // Final discontinuity sensor is just maximum
  DG_FP finalDisSen = fmax(fmax(disSenU_L, disSenV_L), fmax(disSenU_R, disSenV_R));

  out[0][0] += finalDisSen;
  out[1][0] += finalDisSen;

  out_count[0][0] += 1;
  out_count[1][0] += 1;
}
