//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf2_faces_vis_openacc( const double *nuL, const double *rhoL, const double *uL,
                                  const double *opL, double *rhsL,
                                  const double *nuR, const double *rhoR, const double *uR,
                                  const double *opR, double *rhsR) {
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opL[ind + n] * uR[n];
    }
    rhsL[m] += nuL[m] * val / rhoL[m];
  }

  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opR[ind + n] * uL[n];
    }
    rhsR[m] += nuR[m] * val / rhoR[m];
  }
}

// host stub function
void op_par_loop_poisson_mf2_faces_vis(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(33);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[33].name      = name;
  OP_kernels[33].count    += 1;

  int  ninds   = 4;
  int  inds[10] = {0,1,2,-1,3,0,1,2,-1,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf2_faces_vis\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_33
    int part_size = OP_PART_SIZE_33;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map0 = arg0.map_data_d;

    double* data3 = (double*)arg3.data_d;
    double* data8 = (double*)arg8.data_d;
    double *data0 = (double *)arg0.data_d;
    double *data1 = (double *)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data4 = (double *)arg4.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map0,data3,data8,data0,data1,data2,data4)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map0idx;
        int map5idx;
        map0idx = map0[n + set_size1 * 0];
        map5idx = map0[n + set_size1 * 1];


        poisson_mf2_faces_vis_openacc(
          &data0[15 * map0idx],
          &data1[15 * map0idx],
          &data2[15 * map0idx],
          &data3[225 * n],
          &data4[15 * map0idx],
          &data0[15 * map5idx],
          &data1[15 * map5idx],
          &data2[15 * map5idx],
          &data8[225 * n],
          &data4[15 * map5idx]);
      }

    }
    OP_kernels[33].transfer  += Plan->transfer;
    OP_kernels[33].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[33].time     += wall_t2 - wall_t1;
}
