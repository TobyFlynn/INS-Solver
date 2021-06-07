//
// auto-generated by op2.py
//

void ls_advec_edges_omp4_kernel(
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size]) \
    map(to: FMASK_ompkernel[:15])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data3[0:dat3size],data5[0:dat5size],data7[0:dat7size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    int map2idx;
    map1idx = map1[n_op + set_size1 * 0];
    map2idx = map1[n_op + set_size1 * 1];

    const double* arg1_vec[] = {
       &data1[3 * map1idx],
       &data1[3 * map2idx]};
    const double* arg3_vec[] = {
       &data3[3 * map1idx],
       &data3[3 * map2idx]};
    const double* arg5_vec[] = {
       &data5[15 * map1idx],
       &data5[15 * map2idx]};
    double* arg7_vec[] = {
       &data7[15 * map1idx],
       &data7[15 * map2idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const double **x = arg1_vec;
    const double **y = arg3_vec;
    const double **q = arg5_vec;
    double **exQ = arg7_vec;

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse;

    if(edgeR == 0) {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
      } else {
        reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
      }
    } else if(edgeR == 1) {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
      } else {
        reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
      }
    } else {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
      } else {
        reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
      }
    }

    int exInd = 0;
    if(edgeL == 1) exInd = 5;
    else if(edgeL == 2) exInd = 2 * 5;

    int *fmask;

    if(edgeR == 0) {
      fmask = FMASK_ompkernel;
    } else if(edgeR == 1) {
      fmask = &FMASK_ompkernel[5];
    } else {
      fmask = &FMASK_ompkernel[2 * 5];
    }

    for(int i = 0; i < 5; i++) {
      int rInd;
      if(reverse) {
        rInd = fmask[5 - i - 1];
      } else {
        rInd = fmask[i];
      }
      exQ[0][exInd + i] += q[1][rInd];
    }

    exInd = 0;
    if(edgeR == 1) exInd = 5;
    else if(edgeR == 2) exInd = 2 * 5;

    if(edgeL == 0) {
      fmask = FMASK_ompkernel;
    } else if(edgeL == 1) {
      fmask = &FMASK_ompkernel[5];
    } else {
      fmask = &FMASK_ompkernel[2 * 5];
    }

    for(int i = 0; i < 5; i++) {
      int lInd;
      if(reverse) {
        lInd = fmask[5 - i - 1];
      } else {
        lInd = fmask[i];
      }
      exQ[1][exInd + i] += q[0][lInd];
    }
    //end inline func
  }

}
