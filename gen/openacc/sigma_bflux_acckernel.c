//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void sigma_bflux_openacc( const int *bedgeNum, const double *sJ, const double *nx,
                        const double *ny, const double *s, double *sigFx,
                        double *sigFy) {
  int exInd = *bedgeNum * 6;

  for(int i = 0; i < 6; i++) {
    sigFx[exInd + i] += gaussW_g[i] * sJ[exInd + i] * nx[exInd + i] * s[exInd + i];
    sigFy[exInd + i] += gaussW_g[i] * sJ[exInd + i] * ny[exInd + i] * s[exInd + i];
  }
}

// host stub function
void op_par_loop_sigma_bflux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(58);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[58].name      = name;
  OP_kernels[58].count    += 1;

  int  ninds   = 6;
  int  inds[7] = {-1,0,1,2,3,4,5};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_bflux\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_58
    int part_size = OP_PART_SIZE_58;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map1 = arg1.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double *data1 = (double *)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data1,data2,data3,data4,data5,data6)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map1idx;
        map1idx = map1[n + set_size1 * 0];


        sigma_bflux_openacc(
          &data0[1 * n],
          &data1[18 * map1idx],
          &data2[18 * map1idx],
          &data3[18 * map1idx],
          &data4[18 * map1idx],
          &data5[18 * map1idx],
          &data6[18 * map1idx]);
      }

    }
    OP_kernels[58].transfer  += Plan->transfer;
    OP_kernels[58].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[58].time     += wall_t2 - wall_t1;
}
