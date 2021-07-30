//
// auto-generated by op2.py
//

void poisson_cells_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *u = &data0[6*n_op];
    const double *op = &data1[36*n_op];
    double *rhs = &data2[6*n_op];

    //inline function
    
    for(int m = 0; m < 6; m++) {
      int ind = m * 6;
      rhs[m] = 0.0;
      for(int n = 0; n < 6; n++) {
        rhs[m] += op[ind + n] * u[n];
      }
    }
    //end inline func
  }

}
