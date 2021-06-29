//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void pressure_solve_apply_bc_openacc( const int *edgeType, const int *edgeNum,
                             const int *d0, const int *d1, const int *d2,
                             const double *mD0, const double *mD1,
                             const double *mD2, const double *sJ,
                             const double *h, const double *tau, const double *gRho, const double *rho,
                             const double *bc, double *b) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;


  const double *mD, *gVM;
  if(*edgeNum == 0) {
    mD  = mD0;
    gVM = gFInterp0_g;
  } else if(*edgeNum == 1) {
    mD  = mD1;
    gVM = gFInterp1_g;
  } else {
    mD  = mD2;
    gVM = gFInterp2_g;
  }

  double tauA[7];
  for(int i = 0; i < 7; i++) {
    int ind = *edgeNum  * 7 + i;
    tauA[i] = 10 * 0.5 * 5 * 6 * (*h / gRho[ind]);
  }

  double op[7 * 15];


  for(int i = 0; i < 7 * 15; i++) {
    int indT = (i % 7) * 15 + i / 7;
    int indSJ = *edgeNum * 7 + (i % 7);
    int indRho = (i / 7);



    op[i] = gVM[indT] * gaussW_g[i % 7] * sJ[indSJ] * tauA[i % 7]
            - (1.0 / gRho[indSJ]) * mD[indT] * gaussW_g[i % 7] * sJ[indSJ];


  }

  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 7; j++) {
      int op_ind = i * 7 + j;
      int bc_ind = *edgeNum * 7 + j;
      b[i] += op[op_ind] * bc[bc_ind];
    }
  }
}

// host stub function
void op_par_loop_pressure_solve_apply_bc(char const *name, op_set set,
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
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 15;
  op_arg args[15];

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
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(26);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[26].name      = name;
  OP_kernels[26].count    += 1;

  int  ninds   = 10;
  int  inds[15] = {-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_solve_apply_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_26
    int part_size = OP_PART_SIZE_26;
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
    int *map5 = arg5.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data9 = (double *)arg9.data_d;
    double *data10 = (double *)arg10.data_d;
    double *data11 = (double *)arg11.data_d;
    double *data12 = (double *)arg12.data_d;
    double *data13 = (double *)arg13.data_d;
    double *data14 = (double *)arg14.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map5,data0,data1,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map5idx;
        map5idx = map5[n + set_size1 * 0];


        pressure_solve_apply_bc_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &arg3_l,
          &arg4_l,
          &data5[105 * map5idx],
          &data6[105 * map5idx],
          &data7[105 * map5idx],
          &data8[21 * map5idx],
          &data9[1 * map5idx],
          &data10[3 * map5idx],
          &data11[21 * map5idx],
          &data12[15 * map5idx],
          &data13[21 * map5idx],
          &data14[15 * map5idx]);
      }

    }
    OP_kernels[26].transfer  += Plan->transfer;
    OP_kernels[26].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[26].time     += wall_t2 - wall_t1;
}
