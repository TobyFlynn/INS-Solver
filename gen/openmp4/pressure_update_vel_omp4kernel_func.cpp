//
// auto-generated by op2.py
//

void pressure_update_vel_omp4_kernel(
  double *arg0,
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
  double *data11,
  int dat11size,
  int count,
  int num_teams,
  int nthread){

  double arg0_l = *arg0;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *factor = &arg0_l;
    const double *rho = &data1[6*n_op];
    const double *dpdx = &data2[6*n_op];
    const double *dpdy = &data3[6*n_op];
    const double *qt0 = &data4[6*n_op];
    const double *qt1 = &data5[6*n_op];
    double *qtt0 = &data6[6*n_op];
    double *qtt1 = &data7[6*n_op];
    double *dpdn = &data8[9*n_op];
    double *prBC = &data9[12*n_op];
    double *pX = &data10[9*n_op];
    double *pY = &data11[9*n_op];

    //inline function
    
    for(int i = 0; i < 6; i++) {
      qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
      qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];


      dpdn[i] = 0.0;
    }

    for(int i = 0; i < 12; i++) {
      prBC[i] = 0.0;
    }

    for(int i = 0; i < 3 * 3; i++) {
      pX[i] = 0.0;
      pY[i] = 0.0;
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
