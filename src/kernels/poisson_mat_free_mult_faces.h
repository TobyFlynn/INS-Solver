inline void poisson_mat_free_mult_faces(const int **order, const DG_FP *dr,
                                        const DG_FP *ds, const DG_FP *dt,
                                        const DG_FP *mmF0, const DG_FP *mmF1,
                                        const DG_FP *mmF2, const DG_FP *mmF3,
                                        const int *faceNum, const int *fmaskL_corrected,
                                        const int *fmaskR_corrected, const DG_FP *nx,
                                        const DG_FP *ny, const DG_FP *nz,
                                        const DG_FP *fscale, const DG_FP *sJ,
                                        const DG_FP **rx, const DG_FP **sx,
                                        const DG_FP **tx, const DG_FP **ry,
                                        const DG_FP **sy, const DG_FP **ty,
                                        const DG_FP **rz, const DG_FP **sz,
                                        const DG_FP **tz, const DG_FP **in,
                                        DG_FP **out) {
  const int p = order[0][0];
  const DG_FP *dr_mat = &dr[(p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF0_mat = &mmF0[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF1_mat = &mmF1[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF2_mat = &mmF2[(p - 1) * DG_NP * DG_NP];
  const DG_FP *mmF3_mat = &mmF3[(p - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP *mmFL, *mmFR;
  if(faceNum[0] == 0)
    mmFL = mmF0_mat;
  else if(faceNum[0] == 1)
    mmFL = mmF1_mat;
  else if(faceNum[0] == 2)
    mmFL = mmF2_mat;
  else
    mmFL = mmF3_mat;

  if(faceNum[1] == 0)
    mmFR = mmF0_mat;
  else if(faceNum[1] == 1)
    mmFR = mmF1_mat;
  else if(faceNum[1] == 2)
    mmFR = mmF2_mat;
  else
    mmFR = mmF3_mat;

  const int findL = faceNum[0] * dg_npf;
  const int findR = faceNum[1] * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  DG_FP tmp0[DG_NP], tmp1[DG_NP];
  DG_FP tmp2[DG_NP], tmp3[DG_NP];

  for(int i = 0; i < dg_np; i++) {
    tmp1[i] = 0.0;
    tmp3[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      DG_FP tmp_mat_mm_L = sJ[0] * mmFL[ind] * in[0][j];
      out[0][i] += 0.5 * tmp_mat_mm_L * gtau;
      tmp1[i] += tmp_mat_mm_L;

      DG_FP tmp_mat_mm_R = sJ[1] * mmFR[ind] * in[1][j];
      out[1][i] += 0.5 * tmp_mat_mm_R * gtau;
      tmp3[i] += tmp_mat_mm_R;
    }
  }

  for(int i = 0; i < dg_np; i++) {
    tmp0[i] = 0.0;
    tmp2[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      DG_FP tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);
      tmp0[i] += tmp_mat_L * in[0][j];
      out[0][j] += -0.5 * tmp_mat_L * tmp1[i];

      DG_FP tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
      tmp2[i] += tmp_mat_R * in[1][j];
      out[1][j] += -0.5 * tmp_mat_R * tmp3[i];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      DG_FP tmp_mat_mm_L = -0.5 * sJ[0] * mmFL[ind];
      out[0][i] += tmp_mat_mm_L * tmp0[j];

      DG_FP tmp_mat_mm_R = -0.5 * sJ[1] * mmFR[ind];
      out[1][i] += tmp_mat_mm_R * tmp2[j];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    tmp2[i] = 0.0;
    tmp3[i] = 0.0;
    for(int j = 0; j < dg_npf; j++) {
      // int findL = i + fmaskL[j] * dg_np;
      int findL = DG_MAT_IND(i, fmaskL[j], dg_np, dg_np);
      DG_FP tmp_mat_mm_L = sJ[0] * mmFL[findL] * in[1][fmaskR_corrected[j]];
      out[0][i] -= 0.5 * gtau * tmp_mat_mm_L;
      tmp2[i] += tmp_mat_mm_L;

      // int findR = i + fmaskR[j] * dg_np;
      int findR = DG_MAT_IND(i, fmaskR[j], dg_np, dg_np);
      DG_FP tmp_mat_mm_R = sJ[1] * mmFR[findR] * in[0][fmaskL_corrected[j]];
      out[1][i] -= 0.5 * gtau * tmp_mat_mm_R;
      tmp3[i] += tmp_mat_mm_R;
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    tmp0[fmaskL[i]] = 0.0;
    tmp1[fmaskR[i]] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      // int indL = fmaskL[i] + j * dg_np;
      int indL = DG_MAT_IND(fmaskL[i], j, dg_np, dg_np);
      DG_FP tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[indL] + sx[0][0] * ds_mat[indL] + tx[0][0] * dt_mat[indL]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[indL] + sy[0][0] * ds_mat[indL] + ty[0][0] * dt_mat[indL]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[indL] + sz[0][0] * ds_mat[indL] + tz[0][0] * dt_mat[indL]);

      tmp0[fmaskL[i]] += tmp_mat_L * in[0][j];

      // int indR = fmaskR[i] + j * dg_np;
      int indR = DG_MAT_IND(fmaskR[i], j, dg_np, dg_np);
      DG_FP tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[indR] + sx[1][0] * ds_mat[indR] + tx[1][0] * dt_mat[indR]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[indR] + sy[1][0] * ds_mat[indR] + ty[1][0] * dt_mat[indR]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[indR] + sz[1][0] * ds_mat[indR] + tz[1][0] * dt_mat[indR]);

      tmp1[fmaskR[i]] += tmp_mat_R * in[1][j];
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    for(int j = 0; j < dg_npf; j++) {
      // int findL = fmaskL[i] + fmaskL[j] * dg_np;
      int findL = DG_MAT_IND(fmaskL[i], fmaskL[j], dg_np, dg_np);
      out[0][fmaskL[i]] -= 0.5 * sJ[0] * mmFL[findL] * tmp1[fmaskR_corrected[j]];

      // int findR = fmaskR[i] + fmaskR[j] * dg_np;
      int findR = DG_MAT_IND(fmaskR[i], fmaskR[j], dg_np, dg_np);
      out[1][fmaskR[i]] -= 0.5 * sJ[1] * mmFR[findR] * tmp0[fmaskL_corrected[j]];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      DG_FP tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);
      out[0][j] += 0.5 * tmp_mat_L * tmp2[i];

      DG_FP tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
      out[1][j] += 0.5 * tmp_mat_R * tmp3[i];
    }
  }
}
