inline void poisson_mat_free_mult_faces(const int **order, const double *dr,
                                        const double *ds, const double *dt,
                                        const double *mmF0, const double *mmF1,
                                        const double *mmF2, const double *mmF3,
                                        const int *faceNum, const int *fmaskL_corrected,
                                        const int *fmaskR_corrected, const double *nx,
                                        const double *ny, const double *nz,
                                        const double *fscale, const double *sJ,
                                        const double **rx, const double **sx,
                                        const double **tx, const double **ry,
                                        const double **sy, const double **ty,
                                        const double **rz, const double **sz,
                                        const double **tz, const double **in,
                                        double **out) {
  const int p = order[0][0];
  const double *dr_mat = &dr[(p - 1) * DG_NP * DG_NP];
  const double *ds_mat = &ds[(p - 1) * DG_NP * DG_NP];
  const double *dt_mat = &dt[(p - 1) * DG_NP * DG_NP];
  const double *mmF0_mat = &mmF0[(p - 1) * DG_NP * DG_NP];
  const double *mmF1_mat = &mmF1[(p - 1) * DG_NP * DG_NP];
  const double *mmF2_mat = &mmF2[(p - 1) * DG_NP * DG_NP];
  const double *mmF3_mat = &mmF3[(p - 1) * DG_NP * DG_NP];
  const int dg_np  = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const double *mmFL, *mmFR;
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

  const double gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  double tmp0[DG_NP], tmp1[DG_NP];
  double tmp2[DG_NP], tmp3[DG_NP];

  for(int i = 0; i < dg_np; i++) {
    tmp1[i] = 0.0;
    tmp3[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      double tmp_mat_mm_L = sJ[0] * mmFL[ind] * in[0][j];
      out[0][i] += 0.5 * tmp_mat_mm_L * gtau;
      tmp1[i] += tmp_mat_mm_L;

      double tmp_mat_mm_R = sJ[1] * mmFR[ind] * in[1][j];
      out[1][i] += 0.5 * tmp_mat_mm_R * gtau;
      tmp3[i] += tmp_mat_mm_R;
    }
  }

  for(int i = 0; i < dg_np; i++) {
    tmp0[i] = 0.0;
    tmp2[i] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      double tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);
      tmp0[i] += tmp_mat_L * in[0][j];
      out[0][j] += -0.5 * tmp_mat_L * tmp1[i];

      double tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
      tmp2[i] += tmp_mat_R * in[1][j];
      out[1][j] += -0.5 * tmp_mat_R * tmp3[i];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      double tmp_mat_mm_L = -0.5 * sJ[0] * mmFL[ind];
      out[0][i] += tmp_mat_mm_L * tmp0[j];

      double tmp_mat_mm_R = -0.5 * sJ[1] * mmFR[ind];
      out[1][i] += tmp_mat_mm_R * tmp2[j];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    tmp2[i] = 0.0;
    tmp3[i] = 0.0;
    for(int j = 0; j < dg_npf; j++) {
      int findL = i + fmaskL[j] * dg_np;
      double tmp_mat_mm_L = sJ[0] * mmFL[findL] * in[1][fmaskR_corrected[j]];
      out[0][i] -= 0.5 * gtau * tmp_mat_mm_L;
      tmp2[i] += tmp_mat_mm_L;

      int findR = i + fmaskR[j] * dg_np;
      double tmp_mat_mm_R = sJ[1] * mmFR[findR] * in[0][fmaskL_corrected[j]];
      out[1][i] -= 0.5 * gtau * tmp_mat_mm_R;
      tmp3[i] += tmp_mat_mm_R;
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    tmp0[fmaskL[i]] = 0.0;
    tmp1[fmaskR[i]] = 0.0;
    for(int j = 0; j < dg_np; j++) {
      int indL = fmaskL[i] + j * dg_np;
      double tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[indL] + sx[0][0] * ds_mat[indL] + tx[0][0] * dt_mat[indL]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[indL] + sy[0][0] * ds_mat[indL] + ty[0][0] * dt_mat[indL]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[indL] + sz[0][0] * ds_mat[indL] + tz[0][0] * dt_mat[indL]);

      tmp0[fmaskL[i]] += tmp_mat_L * in[0][j];

      int indR = fmaskR[i] + j * dg_np;
      double tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[indR] + sx[1][0] * ds_mat[indR] + tx[1][0] * dt_mat[indR]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[indR] + sy[1][0] * ds_mat[indR] + ty[1][0] * dt_mat[indR]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[indR] + sz[1][0] * ds_mat[indR] + tz[1][0] * dt_mat[indR]);

      tmp1[fmaskR[i]] += tmp_mat_R * in[1][j];
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    for(int j = 0; j < dg_npf; j++) {
      int findL = fmaskL[i] + fmaskL[j] * dg_np;
      out[0][fmaskL[i]] -= 0.5 * sJ[0] * mmFL[findL] * tmp1[fmaskR_corrected[j]];

      int findR = fmaskR[i] + fmaskR[j] * dg_np;
      out[1][fmaskR[i]] -= 0.5 * sJ[1] * mmFR[findR] * tmp0[fmaskL_corrected[j]];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int ind = i + j * dg_np;
      double tmp_mat_L = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      tmp_mat_L += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      tmp_mat_L += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);
      out[0][j] += 0.5 * tmp_mat_L * tmp2[i];

      double tmp_mat_R = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      tmp_mat_R += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      tmp_mat_R += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
      out[1][j] += 0.5 * tmp_mat_R * tmp3[i];
    }
  }
}
