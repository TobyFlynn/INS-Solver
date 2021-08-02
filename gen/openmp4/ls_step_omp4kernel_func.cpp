//
// auto-generated by op2.py
//

void ls_step_omp4_kernel(
  double *arg0,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  int count,
  int num_teams,
  int nthread){

  double arg0_l = *arg0;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size]) \
    map(to: nu0_ompkernel, nu1_ompkernel, rho0_ompkernel, rho1_ompkernel)
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *alpha = &arg0_l;
    const double *s = &data1[10*n_op];
    double *step = &data2[10*n_op];
    double *nu = &data3[10*n_op];
    double *rho = &data4[10*n_op];

    //inline function
    
    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 10; i++) {
      step[i] = tanh(PI * s[i] / *alpha);
      nu[i] = 0.5 * nu0_ompkernel * (1.0 + step[i]) + 0.5 * nu1_ompkernel * (1.0 - step[i]);
      rho[i] = 0.5 * rho0_ompkernel * (1.0 + step[i]) + 0.5 * rho1_ompkernel * (1.0 - step[i]);
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
