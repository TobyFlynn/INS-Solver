//
// auto-generated by op2.py
//

void poisson_op1_omp4_kernel(
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
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size]) \
    map(to: cubW_g_ompkernel[:36])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *J = &data0[36*n_op];
    const double *Dx = &data1[360*n_op];
    const double *Dy = &data2[360*n_op];
    const double *factor = &data3[36*n_op];
    double *op = &data4[100*n_op];

    //inline function
    
    double tmpX[36 * 10];
    double tmpY[36 * 10];

    for(int m = 0; m < 36; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        tmpX[ind] = J[m] * cubW_g_ompkernel[m] * Dx[ind] * factor[m];
        tmpY[ind] = J[m] * cubW_g_ompkernel[m] * Dy[ind] * factor[m];
      }
    }

    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        op[c_ind] = 0.0;
        for(int k = 0; k < 36; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;
          op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
        }
      }
    }
    //end inline func
  }

}
