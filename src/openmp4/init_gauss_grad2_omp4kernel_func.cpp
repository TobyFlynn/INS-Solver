//
// auto-generated by op2.py
//

void init_gauss_grad2_omp4_kernel(
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
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *nx = &data0[21*n_op];
    const double *ny = &data1[21*n_op];
    const double *Dx0 = &data2[105*n_op];
    const double *Dy0 = &data3[105*n_op];
    const double *Dx1 = &data4[105*n_op];
    const double *Dy1 = &data5[105*n_op];
    const double *Dx2 = &data6[105*n_op];
    const double *Dy2 = &data7[105*n_op];
    double *d0 = &data8[105*n_op];
    double *d1 = &data9[105*n_op];
    double *d2 = &data10[105*n_op];

    //inline function
    
    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        d0[ind] = nx[m] * Dx0[ind] + ny[m] * Dy0[ind];
        d1[ind] = nx[m + 7] * Dx1[ind] + ny[m + 7] * Dy1[ind];
        d2[ind] = nx[m + 14] * Dx2[ind] + ny[m + 14] * Dy2[ind];
      }
    }
    //end inline func
  }

}