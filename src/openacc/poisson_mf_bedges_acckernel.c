//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf_bedges_openacc( const int *bedgeType, const int *bedgeNum,
                              const int *d0, const int *d1, const int *d2,
                              const double *sJ, const double *nx,
                              const double *ny, const double *tau,
                              const double *u, const double *dudx,
                              const double *dudy, double *fluxX, double *fluxY,
                              double *flux) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedgeType == *d0 || *bedgeType == *d1 || *bedgeType == *d2) {
    for(int i = 0; i < 7; i++) {
      flux[exInd + i] += gaussW_g[i] * sJ[exInd + i] * (nx[exInd + i] * dudx[exInd + i] + ny[exInd + i] * dudy[exInd + i] - tau[*bedgeNum] * (u[exInd + i]));
    }
  } else {
    for(int i = 0; i < 7; i++) {
      fluxX[exInd + i] += nx[exInd + i] * gaussW_g[i] * sJ[exInd + i] * u[exInd + i];
      fluxY[exInd + i] += ny[exInd + i] * gaussW_g[i] * sJ[exInd + i] * u[exInd + i];
    }
  }
}

// host stub function
void op_par_loop_poisson_mf_bedges(char const *name, op_set set,
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
  op_timing_realloc(34);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;

  int  ninds   = 10;
  int  inds[15] = {-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf_bedges\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_34
    int part_size = OP_PART_SIZE_34;
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


        poisson_mf_bedges_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &arg3_l,
          &arg4_l,
          &data5[21 * map5idx],
          &data6[21 * map5idx],
          &data7[21 * map5idx],
          &data8[3 * map5idx],
          &data9[21 * map5idx],
          &data10[21 * map5idx],
          &data11[21 * map5idx],
          &data12[21 * map5idx],
          &data13[21 * map5idx],
          &data14[21 * map5idx]);
      }

    }
    OP_kernels[34].transfer  += Plan->transfer;
    OP_kernels[34].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].time     += wall_t2 - wall_t1;
}
