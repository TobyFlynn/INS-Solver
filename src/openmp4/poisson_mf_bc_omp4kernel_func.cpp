//
// auto-generated by op2.py
//

void poisson_mf_bc_omp4_kernel(
  int *data0,
  int dat0size,
  int *data1,
  int dat1size,
  int *arg2,
  int *arg3,
  int *arg4,
  int *map5,
  int map5size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  double *data11,
  int dat11size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  int arg2_l = *arg2;
  int arg3_l = *arg3;
  int arg4_l = *arg4;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: gaussW_g_ompkernel[:7], gFInterp0_g_ompkernel[:105], gFInterp1_g_ompkernel[:105], gFInterp2_g_ompkernel[:105])\
    map(to:col_reord[0:set_size1],map5[0:map5size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map5idx;
    map5idx = map5[n_op + set_size1 * 0];

    //variable mapping
    const int *bedgeType = &data0[1*n_op];
    const int *bedgeNum = &data1[1*n_op];
    const int *d0 = &arg2_l;
    const int *d1 = &arg3_l;
    const int *d2 = &arg4_l;
    const double *mD0 = &data5[105 * map5idx];
    const double *mD1 = &data6[105 * map5idx];
    const double *mD2 = &data7[105 * map5idx];
    const double *sJ = &data8[21 * map5idx];
    const double *tau = &data9[3 * map5idx];
    const double *bc = &data10[21 * map5idx];
    double *rhs = &data11[15 * map5idx];

    //inline function
    
    int exInd = 0;
    if(*bedgeNum == 1) exInd = 7;
    else if(*bedgeNum == 2) exInd = 2 * 7;

    if(*bedgeType == *d0 || *bedgeType == *d1 || *bedgeType == *d2) {

      double op[15 * 7];
      for(int j = 0; j < 7 * 15; j++) {
        int indT = (j % 7) * 15 + (j / 7);
        int col  = j % 7;
        int row  = j / 7;
        double val = gaussW_g_ompkernel[j % 7] * sJ[*bedgeNum * 7 + (j % 7)] * tau[*bedgeNum];
        double mD;
        if(*bedgeNum == 0) {
          val *= gFInterp0_g_ompkernel[indT];
          mD = mD0[indT];
        } else if(*bedgeNum == 1) {
          val *= gFInterp1_g_ompkernel[indT];
          mD = mD1[indT];
        } else {
          val *= gFInterp2_g_ompkernel[indT];
          mD = mD2[indT];
        }
        val -= mD * gaussW_g_ompkernel[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
        op[row * 7 + col] += val;
      }

      for(int m = 0; m < 15; m++) {
        int ind = m * 7;
        double val = 0.0;
        for(int n = 0; n < 7; n++) {
          val += op[ind + n] * bc[exInd + n];
        }
        rhs[m] += val;
      }
    } else {
      double op[15 * 7];
      for(int j = 0; j < 7 * 15; j++) {
        int indT = (j % 7) * 15 + (j / 7);
        int col  = j % 7;
        int row  = j / 7;
        double val = gaussW_g_ompkernel[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
        if(*bedgeNum == 0) {
          val *= gFInterp0_g_ompkernel[indT];
        } else if(*bedgeNum == 1) {
          val *= gFInterp1_g_ompkernel[indT];
        } else {
          val *= gFInterp2_g_ompkernel[indT];
        }
        op[row * 7 + col] += val;
      }

      for(int m = 0; m < 15; m++) {
        int ind = m * 7;
        double val = 0.0;
        for(int n = 0; n < 7; n++) {
          val += op[ind + n] * bc[exInd + n];
        }
        rhs[m] += val;
      }
    }
    //end inline func
  }

  *arg2 = arg2_l;
  *arg3 = arg3_l;
  *arg4 = arg4_l;
}
