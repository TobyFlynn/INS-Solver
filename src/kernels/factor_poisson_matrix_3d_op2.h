inline void factor_poisson_matrix_3d_op2(const int **order, const double *dr,
                                         const double *ds, const double *dt,
                                         const double *mmF0, const double *mmF1,
                                         const double *mmF2, const double *mmF3,
                                         const int *faceNum,
                                         const int *fmaskL_corrected,
                                         const int *fmaskR_corrected,
                                         const double *nx, const double *ny,
                                         const double *nz, const double *fscale,
                                         const double *sJ, const double **rx,
                                         const double **sx, const double **tx,
                                         const double **ry, const double **sy,
                                         const double **ty, const double **rz,
                                         const double **sz, const double **tz,
                                         const double **factor, double *op1L,
                                         double *op1R, double *op2L,
                                         double *op2R) {
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

  double DL[DG_NP * DG_NP], DR[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int ind  = i + j * dg_np;

      DL[ind] = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      DL[ind] += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      DL[ind] += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);
      DL[ind] *= factor[0][i];

      DR[ind] = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      DR[ind] += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      DR[ind] += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
      DR[ind] *= factor[1][i];
    }
  }

  double gtau = 0.0;
  for(int i = 0; i < dg_npf; i++) {
    double tmp = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0] * factor[0][fmaskL[i]], fscale[1] * factor[1][fmaskR_corrected[i]]);
    gtau = fmax(gtau, tmp);
  }
  // Do left face
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = i + j * dg_np;
      double tmp = 0.0;
      for(int k = 0; k < dg_np; k++) {
        int a_ind0 = i + k * dg_np;
        int a_ind1 = i * dg_np + k;
        int b_ind  = j * dg_np + k;
        tmp += -sJ[0] * mmFL[a_ind0] * DL[b_ind] - DL[a_ind1] * sJ[0] * mmFL[b_ind];
      }
      op1L[op_ind] += 0.5 * (gtau * sJ[0] * mmFL[op_ind] + tmp);
    }
  }

  for(int i = 0; i < DG_NP * DG_NP; i++) {
    op2L[i] = 0.0;
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_npf; j++) {
      int op_ind = i + fmaskR_corrected[j] * dg_np;
      int find = i + fmaskL[j] * dg_np;
      op2L[op_ind] -= 0.5 * gtau * sJ[0] * mmFL[find];
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = fmaskL[i] + j * dg_np;
      for(int k = 0; k < dg_npf; k++) {
        int a_ind = fmaskL[i] + fmaskL[k] * dg_np;
        int b_ind = j * dg_np + fmaskR_corrected[k];
        op2L[op_ind] -= 0.5 * sJ[0] * mmFL[a_ind] * -DR[b_ind];
      }
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_npf; j++) {
      int op_ind = i + fmaskR_corrected[j] * dg_np;
      for(int k = 0; k < dg_np; k++) {
        int a_ind = i * dg_np + k;
        int b_ind = fmaskL[j] * dg_np + k;
        op2L[op_ind] -= -0.5 * DL[a_ind] * sJ[0] * mmFL[b_ind];
      }
    }
  }

  // Do right face
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = i + j * dg_np;
      double tmp = 0.0;
      for(int k = 0; k < dg_np; k++) {
        int a_ind0 = i + k * dg_np;
        int a_ind1 = i * dg_np + k;
        int b_ind  = j * dg_np + k;
        tmp += -sJ[1] * mmFR[a_ind0] * DR[b_ind] - DR[a_ind1] * sJ[1] * mmFR[b_ind];
      }
      op1R[op_ind] += 0.5 * (gtau * sJ[1] * mmFR[op_ind] + tmp);
    }
  }

  for(int i = 0; i < DG_NP * DG_NP; i++) {
    op2R[i] = 0.0;
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_npf; j++) {
      int op_ind = i + fmaskL_corrected[j] * dg_np;
      int find = i + fmaskR[j] * dg_np;
      op2R[op_ind] -= 0.5 * gtau * sJ[1] * mmFR[find];
    }
  }

  for(int i = 0; i < dg_npf; i++) {
    for(int j = 0; j < dg_np; j++) {
      int op_ind = fmaskR[i] + j * dg_np;
      for(int k = 0; k < dg_npf; k++) {
        int a_ind = fmaskR[i] + fmaskR[k] * dg_np;
        int b_ind = j * dg_np + fmaskL_corrected[k];
        op2R[op_ind] -= 0.5 * sJ[1] * mmFR[a_ind] * -DL[b_ind];
      }
    }
  }

  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_npf; j++) {
      int op_ind = i + fmaskL_corrected[j] * dg_np;
      for(int k = 0; k < dg_np; k++) {
        int a_ind = i * dg_np + k;
        int b_ind = fmaskR[j] * dg_np + k;
        op2R[op_ind] -= -0.5 * DR[a_ind] * sJ[1] * mmFR[b_ind];
      }
    }
  }
}
