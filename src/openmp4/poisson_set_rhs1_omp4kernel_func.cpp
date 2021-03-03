//
// auto-generated by op2.py
//

void poisson_set_rhs1_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *nx = &data0[15*n_op];
    const double *ny = &data1[15*n_op];
    const double *sJ = &data2[15*n_op];
    const double *bc = &data3[15*n_op];
    double *bcNX = &data4[15*n_op];
    double *bcNY = &data5[15*n_op];
    double *tau = &data6[15*n_op];

    //inline function
    
    for(int i = 0; i < 15; i++) {
      bcNX[i] = sJ[i] * nx[i] * bc[i];
      bcNY[i] = sJ[i] * ny[i] * bc[i];
      tau[i]  = sJ[i] * tau[i];
    }
    //end inline func
  }

}
