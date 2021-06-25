//
// auto-generated by op2.py
//

void pressure_solve_apply_bc_omp4_kernel(
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
  double *data12,
  int dat12size,
  double *data13,
  int dat13size,
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
    map(to:col_reord[0:set_size1],map5[0:map5size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size],data12[0:dat12size],data13[0:dat13size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map5idx;
    map5idx = map5[n_op + set_size1 * 0];

    //variable mapping
    const int *edgeType = &data0[1*n_op];
    const int *edgeNum = &data1[1*n_op];
    const int *d0 = &arg2_l;
    const int *d1 = &arg3_l;
    const int *d2 = &arg4_l;
    const double *mD0 = &data5[105 * map5idx];
    const double *mD1 = &data6[105 * map5idx];
    const double *mD2 = &data7[105 * map5idx];
    const double *sJ = &data8[21 * map5idx];
    const double *h = &data9[1 * map5idx];
    const double *tau = &data10[3 * map5idx];
    const double *rho = &data11[21 * map5idx];
    const double *bc = &data12[21 * map5idx];
    double *b = &data13[15 * map5idx];

    //inline function
    
    if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
      return;


    const double *mD, *gVM;
    if(*edgeNum == 0) {
      mD  = mD0;
      gVM = gFInterp0_g_ompkernel;
    } else if(*edgeNum == 1) {
      mD  = mD1;
      gVM = gFInterp1_g_ompkernel;
    } else {
      mD  = mD2;
      gVM = gFInterp2_g_ompkernel;
    }


    double op[7 * 15];


    for(int i = 0; i < 7 * 15; i++) {
      int indT = (i % 7) * 15 + i / 7;
      int indSJ = *edgeNum * 7 + (i % 7);
      int indRho = *edgeNum * 7 + (i / 7);


      op[i] = gVM[indT] * gaussW_g_ompkernel[i % 7] * sJ[indSJ] * tau[*edgeNum]
              - (1.0 / rho[indRho]) * mD[indT] * gaussW_g_ompkernel[i % 7] * sJ[indSJ];
    }

    for(int i = 0; i < 15; i++) {
      for(int j = 0; j < 7; j++) {
        int op_ind = i * 7 + j;
        int bc_ind = *edgeNum * 7 + j;
        b[i] += op[op_ind] * bc[bc_ind];
      }
    }
    //end inline func
  }

  *arg2 = arg2_l;
  *arg3 = arg3_l;
  *arg4 = arg4_l;
}
