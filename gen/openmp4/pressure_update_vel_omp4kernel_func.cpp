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
  int count,
  int num_teams,
  int nthread){

  double arg0_l = *arg0;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *factor = &arg0_l;
    const double *rho = &data1[10*n_op];
    const double *dpdx = &data2[10*n_op];
    const double *dpdy = &data3[10*n_op];
    const double *qt0 = &data4[10*n_op];
    const double *qt1 = &data5[10*n_op];
    double *qtt0 = &data6[10*n_op];
    double *qtt1 = &data7[10*n_op];
    double *dpdn = &data8[12*n_op];

    //inline function
    
    for(int i = 0; i < 10; i++) {
      qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
      qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];


    }

    for(int i = 0; i < 3 * 4; i++) {
      dpdn[i] = 0.0;
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
