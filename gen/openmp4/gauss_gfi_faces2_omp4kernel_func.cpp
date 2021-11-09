//
// auto-generated by op2.py
//

void gauss_gfi_faces2_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size]) \
    map(to: gFInterp0_g_ompkernel[:60], gFInterp1_g_ompkernel[:60], gFInterp2_g_ompkernel[:60])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    double *gVPL = &data2[60*n_op];
    double *gVPR = &data3[60*n_op];

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;

    const double *gFL, *gFR;
    if(edgeL == 0) {
      gFR = gFInterp0_g_ompkernel;
    } else if(edgeL == 1) {
      gFR = gFInterp1_g_ompkernel;
    } else {
      gFR = gFInterp2_g_ompkernel;
    }

    if(edgeR == 0) {
      gFL = gFInterp0_g_ompkernel;
    } else if(edgeR == 1) {
      gFL = gFInterp1_g_ompkernel;
    } else {
      gFL = gFInterp2_g_ompkernel;
    }

    for(int j = 0; j < 10; j++) {
      for(int i = 0; i < 6; i++) {
        int indL, indR;
        if(!reverse) {
          indL = j * 6 + i;
          indR = j * 6 + i;
        } else {
          indL = j * 6 + i;
          indR = j * 6 + (10 - 1 - i);
        }

        gVPL[indL] = gFL[indR];
        gVPR[indR] = gFR[indL];
      }
    }
    //end inline func
  }

}
