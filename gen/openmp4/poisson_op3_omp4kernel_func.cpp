//
// auto-generated by op2.py
//

void poisson_op3_omp4_kernel(
  int *data0,
  int dat0size,
  int *data1,
  int dat1size,
  int *arg2,
  int *arg3,
  int *arg4,
  double *data5,
  int dat5size,
  int *map6,
  int map6size,
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  int arg2_l = *arg2;
  int arg3_l = *arg3;
  int arg4_l = *arg4;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data5[0:dat5size]) \
    map(to: gaussW_g_ompkernel[:6], gFInterp0_g_ompkernel[:60], gFInterp1_g_ompkernel[:60], gFInterp2_g_ompkernel[:60])\
    map(to:col_reord[0:set_size1],map6[0:map6size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map6idx;
    map6idx = map6[n_op + set_size1 * 0];

    //variable mapping
    const int *edgeType = &data0[1*n_op];
    const int *edgeNum = &data1[1*n_op];
    const int *d0 = &arg2_l;
    const int *d1 = &arg3_l;
    const int *d2 = &arg4_l;
    const double *mD = &data5[60*n_op];
    const double *sJ = &data6[18 * map6idx];
    const double *h = &data7[1 * map6idx];
    const double *gFactor = &data8[18 * map6idx];
    const double *factor = &data9[10 * map6idx];
    double *op1 = &data10[100 * map6idx];

    //inline function
    
    if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
      return;

    const double *gVM;
    if(*edgeNum == 0) {
      gVM = gFInterp0_g_ompkernel;
    } else if(*edgeNum == 1) {
      gVM = gFInterp1_g_ompkernel;
    } else {
      gVM = gFInterp2_g_ompkernel;
    }


    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        for(int k = 0; k < 6; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

          int factors_ind = *edgeNum * 6 + k;



          op1[c_ind] += -gVM[a_ind] * gaussW_g_ompkernel[k] * sJ[factors_ind]
                        * gFactor[factors_ind] * mD[b_ind];


        }
      }
    }


    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        for(int k = 0; k < 6; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

          int factors_ind = *edgeNum * 6 + k;



          op1[c_ind] += -gFactor[factors_ind] * mD[a_ind] * gaussW_g_ompkernel[k]
                        * sJ[factors_ind] * gVM[b_ind];




        }
      }
    }

    double tauA[6];
    double maxTau = 0.0;
    for(int i = 0; i < 6; i++) {
      int ind = *edgeNum  * 6 + i;

      tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);


    }






    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        for(int k = 0; k < 6; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

          int factors_ind = *edgeNum * 6 + k;



          op1[c_ind] += gVM[a_ind] * gaussW_g_ompkernel[k] * sJ[factors_ind]
                        * tauA[k] * gVM[b_ind];
        }
      }
    }
    //end inline func
  }

  *arg2 = arg2_l;
  *arg3 = arg3_l;
  *arg4 = arg4_l;
}
