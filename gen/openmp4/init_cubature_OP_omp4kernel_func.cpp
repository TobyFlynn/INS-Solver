//
// auto-generated by op2.py
//

void init_cubature_OP_omp4_kernel(
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
    map(to: cubW_g_ompkernel[:12])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *J = &data0[12*n_op];
    const double *Dx = &data1[36*n_op];
    const double *Dy = &data2[36*n_op];
    double *temp = &data3[36*n_op];
    double *temp2 = &data4[36*n_op];

    //inline function
    
    for(int m = 0; m < 12; m++) {
      for(int n = 0; n < 3; n++) {
        int ind = m * 3 + n;
        temp[ind] = J[m] * cubW_g_ompkernel[m] * Dx[ind];
        temp2[ind] = J[m] * cubW_g_ompkernel[m] * Dy[ind];
      }
    }
    //end inline func
  }

}
