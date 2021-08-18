//
// auto-generated by op2.py
//

void gauss_grad_faces_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data3,
  int dat3size,
  double *data5,
  int dat5size,
  double *data7,
  int dat7size,
  double *data9,
  int dat9size,
  double *data11,
  int dat11size,
  double *data13,
  int dat13size,
  double *data15,
  int dat15size,
  double *data17,
  int dat17size,
  double *data19,
  int dat19size,
  double *data21,
  int dat21size,
  double *data23,
  int dat23size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data3[0:dat3size],data5[0:dat5size],data7[0:dat7size],data9[0:dat9size],data11[0:dat11size],data13[0:dat13size],data15[0:dat15size],data17[0:dat17size],data19[0:dat19size],data21[0:dat21size],data23[0:dat23size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    int map2idx;
    map1idx = map1[n_op + set_size1 * 0];
    map2idx = map1[n_op + set_size1 * 1];

    const double* arg1_vec[] = {
       &data1[60 * map1idx],
       &data1[60 * map2idx]};
    const double* arg3_vec[] = {
       &data3[60 * map1idx],
       &data3[60 * map2idx]};
    const double* arg5_vec[] = {
       &data5[60 * map1idx],
       &data5[60 * map2idx]};
    const double* arg7_vec[] = {
       &data7[60 * map1idx],
       &data7[60 * map2idx]};
    const double* arg9_vec[] = {
       &data9[60 * map1idx],
       &data9[60 * map2idx]};
    const double* arg11_vec[] = {
       &data11[60 * map1idx],
       &data11[60 * map2idx]};
    double* arg13_vec[] = {
       &data13[60 * map1idx],
       &data13[60 * map2idx]};
    double* arg15_vec[] = {
       &data15[60 * map1idx],
       &data15[60 * map2idx]};
    double* arg17_vec[] = {
       &data17[60 * map1idx],
       &data17[60 * map2idx]};
    double* arg19_vec[] = {
       &data19[60 * map1idx],
       &data19[60 * map2idx]};
    double* arg21_vec[] = {
       &data21[60 * map1idx],
       &data21[60 * map2idx]};
    double* arg23_vec[] = {
       &data23[60 * map1idx],
       &data23[60 * map2idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const double **mDx0 = arg1_vec;
    const double **mDy0 = arg3_vec;
    const double **mDx1 = arg5_vec;
    const double **mDy1 = arg7_vec;
    const double **mDx2 = arg9_vec;
    const double **mDy2 = arg11_vec;
    double **pDx0 = arg13_vec;
    double **pDy0 = arg15_vec;
    double **pDx1 = arg17_vec;
    double **pDy1 = arg19_vec;
    double **pDx2 = arg21_vec;
    double **pDy2 = arg23_vec;

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];

    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int indL = m * 10 + n;
        int indR = m * 10 + n;

        if(edgeL == 0) {
          if(edgeR == 0) {
            pDx0[0][indL] += mDx0[1][indR];
            pDy0[0][indL] += mDy0[1][indR];
            pDx0[1][indR] += mDx0[0][indL];
            pDy0[1][indR] += mDy0[0][indL];
          } else if(edgeR == 1) {
            pDx0[0][indL] += mDx1[1][indR];
            pDy0[0][indL] += mDy1[1][indR];
            pDx1[1][indR] += mDx0[0][indL];
            pDy1[1][indR] += mDy0[0][indL];
          } else {
            pDx0[0][indL] += mDx2[1][indR];
            pDy0[0][indL] += mDy2[1][indR];
            pDx2[1][indR] += mDx0[0][indL];
            pDy2[1][indR] += mDy0[0][indL];
          }
        } else if(edgeL == 1) {
          if(edgeR == 0) {
            pDx1[0][indL] += mDx0[1][indR];
            pDy1[0][indL] += mDy0[1][indR];
            pDx0[1][indR] += mDx1[0][indL];
            pDy0[1][indR] += mDy1[0][indL];
          } else if(edgeR == 1) {
            pDx1[0][indL] += mDx1[1][indR];
            pDy1[0][indL] += mDy1[1][indR];
            pDx1[1][indR] += mDx1[0][indL];
            pDy1[1][indR] += mDy1[0][indL];
          } else {
            pDx1[0][indL] += mDx2[1][indR];
            pDy1[0][indL] += mDy2[1][indR];
            pDx2[1][indR] += mDx1[0][indL];
            pDy2[1][indR] += mDy1[0][indL];
          }
        } else {
          if(edgeR == 0) {
            pDx2[0][indL] += mDx0[1][indR];
            pDy2[0][indL] += mDy0[1][indR];
            pDx0[1][indR] += mDx2[0][indL];
            pDy0[1][indR] += mDy2[0][indL];
          } else if(edgeR == 1) {
            pDx2[0][indL] += mDx1[1][indR];
            pDy2[0][indL] += mDy1[1][indR];
            pDx1[1][indR] += mDx2[0][indL];
            pDy1[1][indR] += mDy2[0][indL];
          } else {
            pDx2[0][indL] += mDx2[1][indR];
            pDy2[0][indL] += mDy2[1][indR];
            pDx2[1][indR] += mDx2[0][indL];
            pDy2[1][indR] += mDy2[0][indL];
          }
        }
      }
    }
    //end inline func
  }

}