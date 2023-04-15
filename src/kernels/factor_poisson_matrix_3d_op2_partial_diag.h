inline void factor_poisson_matrix_3d_op2_partial_diag(const int **order,
                  const DG_FP *dr, const DG_FP *ds, const DG_FP *dt,
                  const DG_FP *mmF0, const DG_FP *mmF1, const DG_FP *mmF2,
                  const DG_FP *mmF3, const int *faceNum,
                  const int *fmaskL_corrected, const int *fmaskR_corrected,
                  const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                  const DG_FP *fscale, const DG_FP *sJ, const DG_FP **rx,
                  const DG_FP **sx, const DG_FP **tx, const DG_FP **ry,
                  const DG_FP **sy, const DG_FP **ty, const DG_FP **rz,
                  const DG_FP **sz, const DG_FP **tz, const DG_FP **factor,
                  DG_FP *diagL, DG_FP *diagR) {
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

  const DG_FP tau_order = (DG_FP) p; // (DG_FP) DG_ORDER;

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

  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  const int findL = faceNumL * dg_npf;
  const int findR = faceNumR * dg_npf;
  const int *fmask  = &FMASK[(p - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNumL * dg_npf];
  const int *fmaskR = &fmask[faceNumR * dg_npf];

  DG_FP D[DG_NP * DG_NP];
  const DG_FP r_fact_0 = nx[0] * rx[0][0] + ny[0] * ry[0][0] + nz[0] * rz[0][0];
  const DG_FP s_fact_0 = nx[0] * sx[0][0] + ny[0] * sy[0][0] + nz[0] * sz[0][0];
  const DG_FP t_fact_0 = nx[0] * tx[0][0] + ny[0] * ty[0][0] + nz[0] * tz[0][0];
  for(int j = 0; j < dg_np; j++) {
    #pragma omp simd
    for(int i = 0; i < dg_np; i++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = r_fact_0 * dr_mat[ind] + s_fact_0 * ds_mat[ind] + t_fact_0 * dt_mat[ind];
      D[ind] *= factor[0][i];
    }
  }

  DG_FP gtau = 0.0;
  for(int i = 0; i < dg_npf; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    DG_FP tmp = 2.0 * (tau_order + 1) * (tau_order + 2) * fmax(fscale[0] * factor[0][fmaskL_ind], fscale[1] * factor[1][fmaskR_ind]);
    gtau = fmax(gtau, tmp);
  }

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFL[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFL[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i,  dg_np,  dg_np);
    diagL[i] += 0.5 * sJ[0] * (gtau * mmFL[ind_t3] - tmp);
  }

  const DG_FP r_fact_1 = nx[1] * rx[1][0] + ny[1] * ry[1][0] + nz[1] * rz[1][0];
  const DG_FP s_fact_1 = nx[1] * sx[1][0] + ny[1] * sy[1][0] + nz[1] * sz[1][0];
  const DG_FP t_fact_1 = nx[1] * tx[1][0] + ny[1] * ty[1][0] + nz[1] * tz[1][0];
  for(int j = 0; j < dg_np; j++) {
    #pragma omp simd
    for(int i = 0; i < dg_np; i++) {
      // int ind = i + j * dg_np;
      int ind = DG_MAT_IND(i, j, dg_np, dg_np);
      D[ind] = r_fact_1 * dr_mat[ind] + s_fact_1 * ds_mat[ind] + t_fact_1 * dt_mat[ind];
      D[ind] *= factor[1][i];
    }
  }

  for(int i = 0; i < dg_np; i++) {
    DG_FP tmp = 0.0;
    for(int k = 0; k < dg_np; k++) {
      int a_ind_t1 = DG_MAT_IND(i, k, dg_np, dg_np);
      int b_ind_t1 = DG_MAT_IND(k, i, dg_np, dg_np);
      int a_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      int b_ind_t2 = DG_MAT_IND(k, i, dg_np, dg_np);
      tmp += mmFR[a_ind_t1] * D[b_ind_t1] + D[a_ind_t2] * mmFR[b_ind_t2];
    }
    int ind_t3 = DG_MAT_IND(i, i, dg_np, dg_np);
    diagR[i] += 0.5 * sJ[1] * (gtau * mmFR[ind_t3] - tmp);
  }
}
