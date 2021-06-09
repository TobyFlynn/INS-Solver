//
// auto-generated by op2.py
//

void init_edges_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  bool *data5,
  int dat5size,
  double *data1,
  int dat1size,
  double *data3,
  int dat3size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data5[0:dat5size])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data3[0:dat3size])
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
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const double **x = arg1_vec;
    const double **y = arg3_vec;
    bool *reverse = &data5[1*n_op];

    //inline function
    
    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];

    if(edgeR == 0) {
      if(edgeL == 0) {
        *reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
      } else if(edgeL == 1) {
        *reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
      } else {
        *reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
      }
    } else if(edgeR == 1) {
      if(edgeL == 0) {
        *reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
      } else if(edgeL == 1) {
        *reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
      } else {
        *reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
      }
    } else {
      if(edgeL == 0) {
        *reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
      } else if(edgeL == 1) {
        *reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
      } else {
        *reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
      }
    }
    //end inline func
  }

}