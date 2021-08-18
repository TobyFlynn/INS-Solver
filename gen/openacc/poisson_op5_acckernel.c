//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_op5_openacc( const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD, const double *sJ, const double *h,
                        const double *gFactor, const double *factor,
                        double *op) {

  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g;
  } else {
    gVM = gFInterp2_g;
  }

  for(int i = 0; i < 6 * 10; i++) {
    op[i] = 0.0;
  }

  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2) {


    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      op[i] = gVM[indT] * gaussW_g[i % 6] * sJ[indSJ];
    }
  } else {

    double tauA[6];
    double maxTau = 0.0;
    for(int i = 0; i < 6; i++) {
      int ind = *edgeNum  * 6 + i;

      tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);


    }





    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      int indFactor = (i / 6);

      op[i] = gVM[indT] * gaussW_g[i % 6] * sJ[indSJ] * tauA[i % 6]
              - factor[indFactor] * mD[indT] * gaussW_g[i % 6] * sJ[indSJ];


    }
  }
}

// host stub function
void op_par_loop_poisson_op5(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(22);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[22].name      = name;
  OP_kernels[22].count    += 1;

  int  ninds   = 4;
  int  inds[11] = {-1,-1,-1,-1,-1,-1,0,1,2,3,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op5\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_22
    int part_size = OP_PART_SIZE_22;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  int arg2_l = arg2h[0];
  int arg3_l = arg3h[0];
  int arg4_l = arg4h[0];

  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map6 = arg6.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data10 = (double*)arg10.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data9 = (double *)arg9.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map6,data0,data1,data5,data10,data6,data7,data8,data9)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map6idx;
        map6idx = map6[n + set_size1 * 0];


        poisson_op5_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &arg3_l,
          &arg4_l,
          &data5[60 * n],
          &data6[18 * map6idx],
          &data7[1 * map6idx],
          &data8[18 * map6idx],
          &data9[10 * map6idx],
          &data10[60 * n]);
      }

    }
    OP_kernels[22].transfer  += Plan->transfer;
    OP_kernels[22].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[22].time     += wall_t2 - wall_t1;
}
