inline void ins_3d_shock_cap_art_vis_0(const int *faceNum, const int *fmaskL_corrected,
                    const int *fmaskR_corrected, const DG_FP *nx, const DG_FP *ny, 
                    const DG_FP *nz, const DG_FP *sJ, const DG_FP *fscale, 
                    const DG_FP **h, const DG_FP **u, const DG_FP **v, const DG_FP **w, 
                    DG_FP **out, int **out_count) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;
  const DG_FP *emat = &dg_Emat_kernel[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF * DG_NP];

  DG_FP maxVelU_L = 0.0;
  DG_FP maxVelV_L = 0.0;
  DG_FP maxVelW_L = 0.0;
  DG_FP maxVelU_R = 0.0;
  DG_FP maxVelV_R = 0.0;
  DG_FP maxVelW_R = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    maxVelU_L = fmax(maxVelU_L, fabs(u[0][i]));
    maxVelU_R = fmax(maxVelU_R, fabs(u[1][i]));
    maxVelV_L = fmax(maxVelV_L, fabs(v[0][i]));
    maxVelV_R = fmax(maxVelV_R, fabs(v[1][i]));
    maxVelW_L = fmax(maxVelW_L, fabs(w[0][i]));
    maxVelW_R = fmax(maxVelW_R, fabs(w[1][i]));
  }

  DG_FP w_L[DG_NPF], w_R[DG_NPF];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    w_L[i] = 0.0;
    for(int j = 0; j < DG_NPF; j++) {
      int ind = DG_MAT_IND(fmaskL_ind, fIndL + j, DG_NP, DG_NPF * DG_NUM_FACES);
      w_L[i] += emat[ind];
    }

    const int fmaskR_ind = fmaskR[i];
    w_R[i] = 0.0;
    for(int j = 0; j < DG_NPF; j++) {
      int ind = DG_MAT_IND(fmaskR_ind, fIndR + j, DG_NP, DG_NPF * DG_NUM_FACES);
      w_R[i] += emat[ind];
    }
  }

  DG_FP marker[DG_NPF];
  DG_FP lengthL = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmaskL[i];
    const DG_FP tmp = nx[0] * u[0][fmask_ind] + ny[0] * v[0][fmask_ind] + nz[0] * w[0][fmask_ind];
    marker[i] = tmp < 0.0 ? 1.0 : 0.0;
    lengthL += marker[i] * w_L[i];
  }
  lengthL *= sJ[0];
  

  DG_FP disSenU_L = 0.0;
  DG_FP disSenV_L = 0.0;
  DG_FP disSenW_L = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    if(marker[i] == 1.0) {
      disSenU_L += w_L[i] * (u[0][fmaskL_ind] - u[1][fmaskR_ind]);
      disSenV_L += w_L[i] * (v[0][fmaskL_ind] - v[1][fmaskR_ind]);
      disSenW_L += w_L[i] * (w[0][fmaskL_ind] - w[1][fmaskR_ind]);
    }
  }

  DG_FP lengthR = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmaskR[i];
    const DG_FP tmp = nx[1] * u[1][fmask_ind] + ny[1] * v[1][fmask_ind] + nz[1] * w[1][fmask_ind];
    marker[i] = tmp < 0.0 ? 1.0 : 0.0;
    lengthR += marker[i] * w_R[i];
  }
  lengthR *= sJ[1];

  DG_FP disSenU_R = 0.0;
  DG_FP disSenV_R = 0.0;
  DG_FP disSenW_R = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];

    if(marker[i] == 1.0) {
      disSenU_R += w_R[i] * (u[1][fmaskR_ind] - u[0][fmaskL_ind]);
      disSenV_R += w_R[i] * (v[1][fmaskR_ind] - v[0][fmaskL_ind]);
      disSenW_R += w_R[i] * (w[1][fmaskR_ind] - w[0][fmaskL_ind]);
    }
  }

  DG_FP pow_order = (DG_ORDER + 1.0) / 2.0;
  disSenU_L = lengthL < 1e-10 ? 0.0 : fabs(disSenU_L * sJ[0] / (lengthL * maxVelU_L * pow(h[0][0], pow_order)));
  disSenV_L = lengthL < 1e-10 ? 0.0 : fabs(disSenV_L * sJ[0] / (lengthL * maxVelV_L * pow(h[0][0], pow_order)));
  disSenW_L = lengthL < 1e-10 ? 0.0 : fabs(disSenW_L * sJ[0] / (lengthL * maxVelW_L * pow(h[0][0], pow_order)));
  disSenU_R = lengthR < 1e-10 ? 0.0 : fabs(disSenU_R * sJ[1] / (lengthR * maxVelU_R * pow(h[1][0], pow_order)));
  disSenV_R = lengthR < 1e-10 ? 0.0 : fabs(disSenV_R * sJ[1] / (lengthR * maxVelV_R * pow(h[1][0], pow_order)));
  disSenW_R = lengthR < 1e-10 ? 0.0 : fabs(disSenW_R * sJ[1] / (lengthR * maxVelW_R * pow(h[1][0], pow_order)));

  if(maxVelU_L < 1e-8) disSenU_L = 0.0;
  if(maxVelU_R < 1e-8) disSenU_R = 0.0;
  if(maxVelV_L < 1e-8) disSenV_L = 0.0;
  if(maxVelV_R < 1e-8) disSenV_R = 0.0;
  if(maxVelW_L < 1e-8) disSenW_L = 0.0;
  if(maxVelW_R < 1e-8) disSenW_R = 0.0;

  // Final discontinuity sensor is just maximum
  DG_FP maxL = fmax(fmax(disSenU_L, disSenV_L), disSenW_L);
  DG_FP maxR = fmax(fmax(disSenU_R, disSenV_R), disSenW_R);
  DG_FP finalDisSen = fmax(maxL, maxR);

  out[0][0] += finalDisSen;
  out[1][0] += finalDisSen;
  out[2][0] += finalDisSen;

  out_count[0][0] += 1;
  out_count[1][0] += 1;
  out_count[2][0] += 1;
}
