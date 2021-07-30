//
// auto-generated by op2.py
//

void poisson_pre_omp4_kernel(
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
    const double *in = &data0[3*n_op];
    const double *pre = &data1[9*n_op];
    double *out = &data2[3*n_op];

    //inline function
    
    for(int i = 0; i < 3; i++) {
      out[i] = 0.0;
      for(int j = 0; j < 3; j++) {

        int ind = i * 3 + j;
        out[i] += pre[ind] * in[j];
      }
    }
    //end inline func
  }

}
