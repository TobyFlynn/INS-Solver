//
// auto-generated by op2.py
//

void poisson_op4_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size]) \
    map(to: cubW_g_ompkernel[:36], cubV_g_ompkernel[:360])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *cJ = &data0[36*n_op];
    const double *factor = &data1[36*n_op];
    double *op = &data2[100*n_op];
    double *tmp = &data3[100*n_op];

    //inline function
    
    double cTmp[36 * 10];
    double mm[10 * 10];
    for(int m = 0; m < 36; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        cTmp[ind] = factor[m] * cJ[m] * cubW_g_ompkernel[m] * cubV_g_ompkernel[ind];
      }
    }

    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        mm[c_ind] = 0.0;
        for(int k = 0; k < 36; k++) {
          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

          mm[c_ind] += cubV_g_ompkernel[b_ind] * cTmp[a_ind];
        }
      }
    }

    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        op[c_ind] += mm[c_ind];
        tmp[c_ind] = mm[c_ind];



      }
    }
    //end inline func
  }

}
