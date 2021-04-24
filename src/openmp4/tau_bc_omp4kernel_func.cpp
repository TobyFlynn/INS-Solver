//
// auto-generated by op2.py
//

void tau_bc_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size]) \
    map(to: FMASK_ompkernel[:15])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    map1idx = map1[n_op + set_size1 * 0];

    //variable mapping
    const int *bedgeNum = &data0[1*n_op];
    const double *J = &data1[15 * map1idx];
    const double *sJ = &data2[15 * map1idx];
    double *tau = &data3[3 * map1idx];

    //inline function
    
    double h = 2.0 * J[FMASK_ompkernel[*bedgeNum * 5]] / sJ[*bedgeNum * 5];
    tau[*bedgeNum] += 100.0 * 15.0 / h;
    //end inline func
  }

}
