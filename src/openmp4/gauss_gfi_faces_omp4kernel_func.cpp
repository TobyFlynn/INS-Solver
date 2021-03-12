//
// auto-generated by op2.py
//

void gauss_gfi_faces_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data3,
  int dat3size,
  double *data5,
  int dat5size,
  double *data7,
  int dat7size,
  double *data9,
  int dat9size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size]) \
    map(to: gFInterp0_ompkernel[:105], gFInterp1_ompkernel[:105], gFInterp2_ompkernel[:105])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data3[0:dat3size],data5[0:dat5size],data7[0:dat7size],data9[0:dat9size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    int map2idx;
    map1idx = map1[n_op + set_size1 * 0];
    map2idx = map1[n_op + set_size1 * 1];

    const double* arg1_vec[] = {
       &data1[3 * map1idx],
       &data1[3 * map2idx]};
    const double* arg3_vec[] = {
       &data3[3 * map1idx],
       &data3[3 * map2idx]};
    double* arg5_vec[] = {
       &data5[105 * map1idx],
       &data5[105 * map2idx]};
    double* arg7_vec[] = {
       &data7[105 * map1idx],
       &data7[105 * map2idx]};
    double* arg9_vec[] = {
       &data9[105 * map1idx],
       &data9[105 * map2idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const double **x = arg1_vec;
    const double **y = arg3_vec;
    double **gf0 = arg5_vec;
    double **gf1 = arg7_vec;
    double **gf2 = arg9_vec;

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse;

    if(edgeR == 0) {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
      } else {
        reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
      }
    } else if(edgeR == 1) {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
      } else {
        reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
      }
    } else {
      if(edgeL == 0) {
        reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
      } else if(edgeL == 1) {
        reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
      } else {
        reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
      }
    }

    for(int m = 0; m < 7; m++) {
      for(int n = 0; n < 15; n++) {
        int indL, indR;
        if(!reverse) {
          indL = m * 15 + n;
          indR = m * 15 + n;
        } else {
          indL = m * 15 + n;
          indR = (6 - m) * 15 + n;
        }

        if(edgeL == 0) {
          if(edgeR == 0) {
            gf0[0][indL] += gFInterp0_ompkernel[indR];
            gf0[1][indR] += gFInterp0_ompkernel[indL];
          } else if(edgeR == 1) {
            gf0[0][indL] += gFInterp1_ompkernel[indR];
            gf1[1][indR] += gFInterp0_ompkernel[indL];
          } else {
            gf0[0][indL] += gFInterp2_ompkernel[indR];
            gf2[1][indR] += gFInterp0_ompkernel[indL];
          }
        } else if(edgeL == 1) {
          if(edgeR == 0) {
            gf1[0][indL] += gFInterp0_ompkernel[indR];
            gf0[1][indR] += gFInterp1_ompkernel[indL];
          } else if(edgeR == 1) {
            gf1[0][indL] += gFInterp1_ompkernel[indR];
            gf1[1][indR] += gFInterp1_ompkernel[indL];
          } else {
            gf1[0][indL] += gFInterp2_ompkernel[indR];
            gf2[1][indR] += gFInterp1_ompkernel[indL];
          }
        } else {
          if(edgeR == 0) {
            gf2[0][indL] += gFInterp0_ompkernel[indR];
            gf0[1][indR] += gFInterp2_ompkernel[indL];
          } else if(edgeR == 1) {
            gf2[0][indL] += gFInterp1_ompkernel[indR];
            gf1[1][indR] += gFInterp2_ompkernel[indL];
          } else {
            gf2[0][indL] += gFInterp2_ompkernel[indR];
            gf2[1][indR] += gFInterp2_ompkernel[indL];
          }
        }
      }
    }
    //end inline func
  }

}
