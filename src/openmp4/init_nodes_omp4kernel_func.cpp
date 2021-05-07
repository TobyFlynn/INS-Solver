//
// auto-generated by op2.py
//

void init_nodes_omp4_kernel(
  int *map0,
  int map0size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data0,
  int dat0size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data3[0:dat3size],data4[0:dat4size])\
    map(to:col_reord[0:set_size1],map0[0:map0size],data0[0:dat0size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map0idx;
    int map1idx;
    int map2idx;
    map0idx = map0[n_op + set_size1 * 0];
    map1idx = map0[n_op + set_size1 * 1];
    map2idx = map0[n_op + set_size1 * 2];

    const double* arg0_vec[] = {
       &data0[2 * map0idx],
       &data0[2 * map1idx],
       &data0[2 * map2idx]};
    //variable mapping
    const double **nc = arg0_vec;
    double *nodeX = &data3[3*n_op];
    double *nodeY = &data4[3*n_op];

    //inline function
    
    nodeX[0] = nc[0][0];
    nodeX[1] = nc[1][0];
    nodeX[2] = nc[2][0];
    nodeY[0] = nc[0][1];
    nodeY[1] = nc[1][1];
    nodeY[2] = nc[2][1];
    //end inline func
  }

}