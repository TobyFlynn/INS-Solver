//
// auto-generated by op2.py
//

void update_Q_omp4_kernel(
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
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *dt = &arg0_l;
    double *s = &data1[3*n_op];
    const double *rk0 = &data2[3*n_op];
    const double *rk1 = &data3[3*n_op];
    const double *rk2 = &data4[3*n_op];

    //inline function
    
    for(int i = 0; i < 3; i++) {
      double add = (*dt) * (rk0[i]/ 6.0 + rk1[i] / 6.0 + 2.0 * rk2[i] / 3.0);
      s[i] = s[i] + add;
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
