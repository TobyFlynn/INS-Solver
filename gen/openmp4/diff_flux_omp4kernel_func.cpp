//
// auto-generated by op2.py
//

void diff_flux_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  int *map2,
  int map2size,
  double *data2,
  int dat2size,
  double *data4,
  int dat4size,
  double *data6,
  int dat6size,
  double *data8,
  int dat8size,
  double *data10,
  int dat10size,
  double *data12,
  int dat12size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: gaussW_g_ompkernel[:3])\
    map(to:col_reord[0:set_size1],map2[0:map2size],data2[0:dat2size],data4[0:dat4size],data6[0:dat6size],data8[0:dat8size],data10[0:dat10size],data12[0:dat12size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map2idx;
    int map3idx;
    map2idx = map2[n_op + set_size1 * 0];
    map3idx = map2[n_op + set_size1 * 1];

    const double* arg2_vec[] = {
       &data2[9 * map2idx],
       &data2[9 * map3idx]};
    const double* arg4_vec[] = {
       &data4[9 * map2idx],
       &data4[9 * map3idx]};
    const double* arg6_vec[] = {
       &data6[9 * map2idx],
       &data6[9 * map3idx]};
    const double* arg8_vec[] = {
       &data8[9 * map2idx],
       &data8[9 * map3idx]};
    const double* arg10_vec[] = {
       &data10[9 * map2idx],
       &data10[9 * map3idx]};
    double* arg12_vec[] = {
       &data12[9 * map2idx],
       &data12[9 * map3idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    const double **sJ = arg2_vec;
    const double **nx = arg4_vec;
    const double **ny = arg6_vec;
    const double **sigX = arg8_vec;
    const double **sigY = arg10_vec;
    double **flux = arg12_vec;

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;

    int exIndL = edgeL * 3;
    int exIndR = edgeR * 3;

    for(int i = 0; i < 3; i++) {
      int rInd;
      int lInd = exIndL + i;
      if(reverse) {
        rInd = exIndR + 3 - i - 1;
      } else {
        rInd = exIndR + i;
      }

      double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
      double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;

      flux[0][lInd] += gaussW_g_ompkernel[i] * sJ[0][lInd] * (nx[0][lInd] * sigFX + ny[0][lInd] * sigFY);
    }

    for(int i = 0; i < 3; i++) {
      int lInd;
      int rInd = exIndR + i;
      if(reverse) {
        lInd = exIndL + 3 - i - 1;
      } else {
        lInd = exIndL + i;
      }

      double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
      double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;

      flux[1][rInd] += gaussW_g_ompkernel[i] * sJ[1][rInd] * (nx[1][rInd] * sigFX + ny[1][rInd] * sigFY);
    }
    //end inline func
  }

}
