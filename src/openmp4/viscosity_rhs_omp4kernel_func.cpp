//
// auto-generated by op2.py
//

void viscosity_rhs_omp4_kernel(
  double *arg0,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  int count,
  int num_teams,
  int nthread){

  double arg0_l = *arg0;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *factor = &arg0_l;
    double *vRHS0 = &data1[15*n_op];
    double *vRHS1 = &data2[15*n_op];

    //inline function
    
    for(int i = 0; i < 15; i++) {
      vRHS0[i] = (*factor) * vRHS0[i];
      vRHS1[i] = (*factor) * vRHS1[i];
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
