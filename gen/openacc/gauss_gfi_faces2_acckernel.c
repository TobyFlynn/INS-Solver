//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void gauss_gfi_faces2_openacc( const int *edgeNum, const bool *rev,
                             double *gVPL, double *gVPR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  const double *gFL, *gFR;
  if(edgeL == 0) {
    gFR = gFInterp0_g;
  } else if(edgeL == 1) {
    gFR = gFInterp1_g;
  } else {
    gFR = gFInterp2_g;
  }

  if(edgeR == 0) {
    gFL = gFInterp0_g;
  } else if(edgeR == 1) {
    gFL = gFInterp1_g;
  } else {
    gFL = gFInterp2_g;
  }

  for(int i = 0; i < 6 * 10; i++) {
    gVPL[i] = 0.0;
    gVPR[i] = 0.0;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 10 + n;
        indR = m * 10 + n;
      } else {
        indL = m * 10 + n;
        indR = (6 - 1 - m) * 10 + n;
      }

      gVPL[indL] += gFL[indR];
      gVPR[indR] += gFR[indL];
    }
  }
}

// host stub function
void op_par_loop_gauss_gfi_faces2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  gauss_gfi_faces2");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    int* data0 = (int*)arg0.data_d;
    bool* data1 = (bool*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3)
    for ( int n=0; n<set->size; n++ ){
      gauss_gfi_faces2_openacc(
        &data0[2*n],
        &data1[1*n],
        &data2[60*n],
        &data3[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg0.size;
  OP_kernels[12].transfer += (float)set->size * arg1.size;
  OP_kernels[12].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[12].transfer += (float)set->size * arg3.size * 2.0f;
}