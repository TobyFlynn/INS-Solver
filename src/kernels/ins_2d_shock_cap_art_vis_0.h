inline void ins_2d_shock_cap_art_vis_0(const int *faceNum, const bool *reverse,
                              const DG_FP *nx, const DG_FP *ny, const DG_FP *sJ,
                              const DG_FP *fscale, const DG_FP **h,
                              const DG_FP **node_coords, const DG_FP **u,
                              const DG_FP **v, DG_FP **out, DG_FP **out_count) {
  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int *fmask  = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * DG_NPF];
  const int *fmaskR = &fmask[faceNumR * DG_NPF];
  const bool rev = *reverse;

  const DG_FP *emat = &dg_Emat_kernel[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF * DG_NP];

  DG_FP maxVelU_L = 0.0;
  DG_FP maxVelV_L = 0.0;
  DG_FP maxVelU_R = 0.0;
  DG_FP maxVelV_R = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    maxVelU_L = fmax(maxVelU_L, fabs(u[0][i]));
    maxVelU_R = fmax(maxVelU_R, fabs(u[1][i]));
    maxVelV_L = fmax(maxVelV_L, fabs(v[0][i]));
    maxVelV_R = fmax(maxVelV_R, fabs(v[1][i]));
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

  DG_FP marker[DG_NPF];
  DG_FP lengthL = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmaskL[i];
    const DG_FP tmp = nx[0] * u[0][fmask_ind] + ny[0] * v[0][fmask_ind];
    marker[i] = tmp < 0.0 ? 1.0 : 0.0;
    lengthL += marker[i] * w_L[i];
  }
  lengthL *= sJ[0];


  DG_FP disSenU_L = 0.0;
  DG_FP disSenV_L = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = rev ? fmaskR[DG_NPF - i - 1] : fmaskR[i];
    if(marker[i] == 1.0) {
      disSenU_L += w_L[i] * (u[0][fmaskL_ind] - u[1][fmaskR_ind]);
      disSenV_L += w_L[i] * (v[0][fmaskL_ind] - v[1][fmaskR_ind]);
    }
  }

  DG_FP lengthR = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmaskR[i];
    const DG_FP tmp = nx[1] * u[1][fmask_ind] + ny[1] * v[1][fmask_ind];
    marker[i] = tmp < 0.0 ? 1.0 : 0.0;
    lengthR += marker[i] * w_R[i];
  }
  lengthR *= sJ[1];

  DG_FP disSenU_R = 0.0;
  DG_FP disSenV_R = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = rev ? fmaskL[DG_NPF - i - 1] : fmaskL[i];

    if(marker[i] == 1.0) {
      disSenU_R += w_R[i] * (u[1][fmaskR_ind] - u[0][fmaskL_ind]);
      disSenV_R += w_R[i] * (v[1][fmaskR_ind] - v[0][fmaskL_ind]);
    }
  }

  // DG_FP x_diff = node_coords[0][0] - node_coords[1][0];
  // DG_FP y_diff = node_coords[0][1] - node_coords[1][1];
  // DG_FP edge_size = sqrt(x_diff * x_diff + y_diff * y_diff);
  DG_FP pow_order = (DG_ORDER + 1.0) / 2.0;
  disSenU_L = lengthL < 1e-10 ? 0.0 : fabs(disSenU_L * sJ[0] / (lengthL * maxVelU_L * pow(h[0][0], pow_order)));
  disSenV_L = lengthL < 1e-10 ? 0.0 : fabs(disSenV_L * sJ[0] / (lengthL * maxVelV_L * pow(h[0][0], pow_order)));
  disSenU_R = lengthR < 1e-10 ? 0.0 : fabs(disSenU_R * sJ[1] / (lengthR * maxVelU_R * pow(h[1][0], pow_order)));
  disSenV_R = lengthR < 1e-10 ? 0.0 : fabs(disSenV_R * sJ[1] / (lengthR * maxVelV_R * pow(h[1][0], pow_order)));
  // disSenU_L = fabs(disSenU_L * sJ[0] / (edge_size * maxVelU_L));
  // disSenV_L = fabs(disSenV_L * sJ[0] / (edge_size * maxVelV_L));
  // disSenU_R = fabs(disSenU_R * sJ[1] / (edge_size * maxVelU_R));
  // disSenV_R = fabs(disSenV_R * sJ[1] / (edge_size * maxVelV_R));

  if(maxVelU_L < 1e-8) disSenU_L = 0.0;
  if(maxVelU_R < 1e-8) disSenU_R = 0.0;
  if(maxVelV_L < 1e-8) disSenV_L = 0.0;
  if(maxVelV_R < 1e-8) disSenV_R = 0.0;

  // Final discontinuity sensor is just maximum
  DG_FP finalDisSen = fmax(fmax(disSenU_L, disSenV_L), fmax(disSenU_R, disSenV_R));

  out[0][0] += finalDisSen;
  out[1][0] += finalDisSen;

  out_count[0][0] += 1.0;
  out_count[1][0] += 1.0;
}
