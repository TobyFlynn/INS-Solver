inline void poisson_matrix_3d_op2_partial(const int **order, const DG_FP *dr,
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
                                  const DG_FP **tz, DG_FP *op1L, DG_FP *op1R) {
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

  DG_FP DL[DG_NP * DG_NP], DR[DG_NP * DG_NP];
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);

      DL[ind] = nx[0] * (rx[0][0] * dr_mat[ind] + sx[0][0] * ds_mat[ind] + tx[0][0] * dt_mat[ind]);
      DL[ind] += ny[0] * (ry[0][0] * dr_mat[ind] + sy[0][0] * ds_mat[ind] + ty[0][0] * dt_mat[ind]);
      DL[ind] += nz[0] * (rz[0][0] * dr_mat[ind] + sz[0][0] * ds_mat[ind] + tz[0][0] * dt_mat[ind]);

      DR[ind] = nx[1] * (rx[1][0] * dr_mat[ind] + sx[1][0] * ds_mat[ind] + tx[1][0] * dt_mat[ind]);
      DR[ind] += ny[1] * (ry[1][0] * dr_mat[ind] + sy[1][0] * ds_mat[ind] + ty[1][0] * dt_mat[ind]);
      DR[ind] += nz[1] * (rz[1][0] * dr_mat[ind] + sz[1][0] * ds_mat[ind] + tz[1][0] * dt_mat[ind]);
    }
  }

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, -0.5 * sJ[0], mmFL, dg_np, DL, dg_np, 1.0, op1L, dg_np);
  op2_in_kernel_gemm(true, false, dg_np, dg_np, dg_np, -0.5 * sJ[0], DL, dg_np, mmFL, dg_np, 1.0, op1L, dg_np);

  const DG_FP tmp_constL = 0.5 * sJ[0] * gtau;
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int op_ind = i + j * dg_np;
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      op1L[op_ind] += tmp_constL * mmFL[op_ind];
    }
  }

  op2_in_kernel_gemm(false, false, dg_np, dg_np, dg_np, -0.5 * sJ[1], mmFR, dg_np, DR, dg_np, 1.0, op1R, dg_np);
  op2_in_kernel_gemm(true, false, dg_np, dg_np, dg_np, -0.5 * sJ[1], DR, dg_np, mmFR, dg_np, 1.0, op1R, dg_np);

  const DG_FP tmp_constR = 0.5 * sJ[1] * gtau;
  for(int i = 0; i < dg_np; i++) {
    for(int j = 0; j < dg_np; j++) {
      // int op_ind = i + j * dg_np;
      int op_ind = DG_MAT_IND(i, j, dg_np, dg_np);
      op1R[op_ind] += tmp_constR * mmFR[op_ind];
    }
  }
}
