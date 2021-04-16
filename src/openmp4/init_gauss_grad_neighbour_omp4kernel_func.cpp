//
// auto-generated by op2.py
//

void init_gauss_grad_neighbour_omp4_kernel(
  int *data0,
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

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size]) \
    map(to: gF0Dr_g_ompkernel[:105], gF0Ds_g_ompkernel[:105], gF1Dr_g_ompkernel[:105], gF1Ds_g_ompkernel[:105], gF2Dr_g_ompkernel[:105], gF2Ds_g_ompkernel[:105], gF0DrR_g_ompkernel[:105], gF0DsR_g_ompkernel[:105], gF1DrR_g_ompkernel[:105], gF1DsR_g_ompkernel[:105], gF2DrR_g_ompkernel[:105], gF2DsR_g_ompkernel[:105])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const int *reverse = &data0[3*n_op];
    double *rx = &data1[21*n_op];
    double *sx = &data2[21*n_op];
    double *ry = &data3[21*n_op];
    double *sy = &data4[21*n_op];
    double *Dx0 = &data5[105*n_op];
    double *Dy0 = &data6[105*n_op];
    double *Dx1 = &data7[105*n_op];
    double *Dy1 = &data8[105*n_op];
    double *Dx2 = &data9[105*n_op];
    double *Dy2 = &data10[105*n_op];

    //inline function
    

    double J[21];
    for(int i = 0; i < 21; i++) {
      J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
    }

    for(int i = 0; i < 21; i++) {
      double rx_n = sy[i] / J[i];
      double sx_n = -ry[i] / J[i];
      double ry_n = -sx[i] / J[i];
      double sy_n = rx[i] / J[i];
      rx[i] = rx_n;
      sx[i] = sx_n;
      ry[i] = ry_n;
      sy[i] = sy_n;
    }

    if(reverse[0]) {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx0[ind] = rx[m] * gF0DrR_g_ompkernel[ind] + sx[m] * gF0DsR_g_ompkernel[ind];
          Dy0[ind] = ry[m] * gF0DrR_g_ompkernel[ind] + sy[m] * gF0DsR_g_ompkernel[ind];
        }
      }
    } else {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx0[ind] = rx[m] * gF0Dr_g_ompkernel[ind] + sx[m] * gF0Ds_g_ompkernel[ind];
          Dy0[ind] = ry[m] * gF0Dr_g_ompkernel[ind] + sy[m] * gF0Ds_g_ompkernel[ind];
        }
      }
    }

    if(reverse[1]) {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx1[ind] = rx[m + 7] * gF1DrR_g_ompkernel[ind] + sx[m + 7] * gF1DsR_g_ompkernel[ind];
          Dy1[ind] = ry[m + 7] * gF1DrR_g_ompkernel[ind] + sy[m + 7] * gF1DsR_g_ompkernel[ind];
        }
      }
    } else {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx1[ind] = rx[m + 7] * gF1Dr_g_ompkernel[ind] + sx[m + 7] * gF1Ds_g_ompkernel[ind];
          Dy1[ind] = ry[m + 7] * gF1Dr_g_ompkernel[ind] + sy[m + 7] * gF1Ds_g_ompkernel[ind];
        }
      }
    }

    if(reverse[2]) {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx2[ind] = rx[m + 14] * gF2DrR_g_ompkernel[ind] + sx[m + 14] * gF2DsR_g_ompkernel[ind];
          Dy2[ind] = ry[m + 14] * gF2DrR_g_ompkernel[ind] + sy[m + 14] * gF2DsR_g_ompkernel[ind];
        }
      }
    } else {
      for(int m = 0; m < 7; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          Dx2[ind] = rx[m + 14] * gF2Dr_g_ompkernel[ind] + sx[m + 14] * gF2Ds_g_ompkernel[ind];
          Dy2[ind] = ry[m + 14] * gF2Dr_g_ompkernel[ind] + sy[m + 14] * gF2Ds_g_ompkernel[ind];
        }
      }
    }
    //end inline func
  }

}
