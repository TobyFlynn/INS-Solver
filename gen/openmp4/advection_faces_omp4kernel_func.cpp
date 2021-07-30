//
// auto-generated by op2.py
//

void advection_faces_omp4_kernel(
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: FMASK_ompkernel[:9])\
    map(to:col_reord[0:set_size1],map2[0:map2size],data2[0:dat2size],data4[0:dat4size],data6[0:dat6size],data8[0:dat8size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map2idx;
    int map3idx;
    map2idx = map2[n_op + set_size1 * 0];
    map3idx = map2[n_op + set_size1 * 1];

    const double* arg2_vec[] = {
       &data2[6 * map2idx],
       &data2[6 * map3idx]};
    const double* arg4_vec[] = {
       &data4[6 * map2idx],
       &data4[6 * map3idx]};
    double* arg6_vec[] = {
       &data6[9 * map2idx],
       &data6[9 * map3idx]};
    double* arg8_vec[] = {
       &data8[9 * map2idx],
       &data8[9 * map3idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    const double **q0 = arg2_vec;
    const double **q1 = arg4_vec;
    double **exQ0 = arg6_vec;
    double **exQ1 = arg8_vec;

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;

    int exInd = edgeL * 3;
    int *fmask = &FMASK_ompkernel[edgeR * 3];

    for(int i = 0; i < 3; i++) {
      int rInd;
      if(reverse) {
        rInd = fmask[3 - i - 1];
      } else {
        rInd = fmask[i];
      }
      exQ0[0][exInd + i] += q0[1][rInd];
      exQ1[0][exInd + i] += q1[1][rInd];
    }

    exInd = edgeR * 3;
    fmask = &FMASK_ompkernel[edgeL * 3];

    for(int i = 0; i < 3; i++) {
      int lInd;
      if(reverse) {
        lInd = fmask[3 - i - 1];
      } else {
        lInd = fmask[i];
      }
      exQ0[1][exInd + i] += q0[0][lInd];
      exQ1[1][exInd + i] += q1[0][lInd];
    }
    //end inline func
  }

}
