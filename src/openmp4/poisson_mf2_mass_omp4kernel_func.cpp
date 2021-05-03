//
// auto-generated by op2.py
//

void poisson_mf2_mass_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *arg2,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  int count,
  int num_teams,
  int nthread){

  double arg2_l = *arg2;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data3[0:dat3size],data4[0:dat4size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *u = &data0[15*n_op];
    const double *op = &data1[225*n_op];
    const double *factor = &arg2_l;
    const double *mm = &data3[225*n_op];
    double *rhs = &data4[15*n_op];

    //inline function
    
    double mFactor = *factor;
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {

        val += (op[ind + n] + mm[n * 15 + m] * mFactor) * u[n];
      }
      rhs[m] = val;
    }
    //end inline func
  }

  *arg2 = arg2_l;
}