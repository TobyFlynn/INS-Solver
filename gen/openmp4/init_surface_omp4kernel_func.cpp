//
// auto-generated by op2.py
//

void init_surface_omp4_kernel(
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
    const double *x = &data0[10*n_op];
    const double *y = &data1[10*n_op];
    double *s = &data2[10*n_op];

    //inline function
    
    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 10; i++) {

      s[i] = sqrt((x[i] - 1.0) * (x[i] - 1.0) + (y[i] - 0.5) * (y[i] - 0.5)) - 0.15;
    }
    //end inline func
  }

}