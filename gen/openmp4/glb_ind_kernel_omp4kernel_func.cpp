//
// auto-generated by op2.py
//

void glb_ind_kernel_omp4_kernel(
  int *map0,
  int map0size,
  int *data2,
  int dat2size,
  int *data3,
  int dat3size,
  int *data0,
  int dat0size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data2[0:dat2size],data3[0:dat3size])\
    map(to:col_reord[0:set_size1],map0[0:map0size],data0[0:dat0size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map0idx;
    int map1idx;
    map0idx = map0[n_op + set_size1 * 0];
    map1idx = map0[n_op + set_size1 * 1];

    const int* arg0_vec[] = {
       &data0[1 * map0idx],
       &data0[1 * map1idx]};
    //variable mapping
    const int **glb = arg0_vec;
    int *glbL = &data2[1*n_op];
    int *glbR = &data3[1*n_op];

    //inline function
    
    glbL[0] = glb[0][0];
    glbR[0] = glb[1][0];
    //end inline func
  }

}