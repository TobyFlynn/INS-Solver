//
// auto-generated by op2.py
//

void poisson_edges_omp4_kernel(
  int *map0,
  int map0size,
  double *data1,
  int dat1size,
  double *data4,
  int dat4size,
  double *data0,
  int dat0size,
  double *data2,
  int dat2size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data4[0:dat4size])\
    map(to:col_reord[0:set_size1],map0[0:map0size],data0[0:dat0size],data2[0:dat2size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map0idx;
    int map3idx;
    map0idx = map0[n_op + set_size1 * 0];
    map3idx = map0[n_op + set_size1 * 1];

    //variable mapping
    const double *uL = &data0[6 * map0idx];
    const double *opL = &data1[36*n_op];
    double *rhsL = &data2[6 * map0idx];
    const double *uR = &data0[6 * map3idx];
    const double *opR = &data4[36*n_op];
    double *rhsR = &data2[6 * map3idx];

    //inline function
    
    for(int m = 0; m < 6; m++) {
      int ind = m * 6;
      for(int n = 0; n < 6; n++) {
        rhsL[m] += opL[ind + n] * uR[n];
        rhsR[m] += opR[ind + n] * uL[n];
      }
    }
    //end inline func
  }

}
