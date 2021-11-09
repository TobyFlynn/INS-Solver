//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void pressure_bc_openacc( const int *bedge_type, const int *bedgeNum,
                        const double *t, const int *problem, const double *x,
                        const double *y, const double *nx, const double *ny,
                        const double *nu, const double *rho, const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        double *dPdN) {
  int exInd = *bedgeNum * 4;
  int *fmask = &FMASK[*bedgeNum * 4];

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0 || *bedge_type == 2 || *bedge_type == 3) {

      for(int i = 0; i < 4; i++) {
        int fInd = fmask[i];


        double res1 = -N0[fInd] - gradCurlVel1[fInd] / (reynolds * rho[fInd]);
        double res2 = -N1[fInd] + gradCurlVel0[fInd] / (reynolds * rho[fInd]);
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
      }
    }

    if(*bedge_type == 0) {

      for(int i = 0; i < 4; i++) {
        double y1 = y[fmask[i]];
        int fInd = fmask[i];
        double bcdUndt = -pow(1.0, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0) * 6.0 * y1 * (1.0 - y1);

        dPdN[exInd + i] -= bcdUndt;
      }
    }
  } else {
    if(*bedge_type == 0) {

      for(int i = 0; i < 4; i++) {
        int fInd = fmask[i];
        double res1 = -N0[fInd] - nu[fInd] * gradCurlVel1[fInd];
        double res2 = -N1[fInd] + nu[fInd] * gradCurlVel0[fInd];
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;

        double y1 = y[fmask[i]];
        double x1 = x[fmask[i]];
        double nx1 = nx[exInd + i];
        double ny1 = ny[exInd + i];
        double bcdUndt = -nu[fInd] * 4.0 * PI * PI * (-nx1 * sin(2.0 * PI * y1) + ny1 * sin(2.0 * PI * x1))
                          * exp(-nu[fInd] * 4.0 * PI * PI * *t);
        dPdN[exInd + i] -= bcdUndt;
      }
    }
  }
}

// host stub function
void op_par_loop_pressure_bc(char const *name, op_set set,
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

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
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
  op_timing_realloc(66);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[66].name      = name;
  OP_kernels[66].count    += 1;

  int  ninds   = 11;
  int  inds[15] = {-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9,10};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_66
    int part_size = OP_PART_SIZE_66;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg2_l = arg2h[0];
  int arg3_l = arg3h[0];

  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map4 = arg4.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data4 = (double *)arg4.data_d;
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

      #pragma acc parallel loop independent deviceptr(col_reord,map4,data0,data1,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map4idx;
        map4idx = map4[n + set_size1 * 0];


        pressure_bc_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &arg3_l,
          &data4[10 * map4idx],
          &data5[10 * map4idx],
          &data6[12 * map4idx],
          &data7[12 * map4idx],
          &data8[10 * map4idx],
          &data9[10 * map4idx],
          &data10[10 * map4idx],
          &data11[10 * map4idx],
          &data12[10 * map4idx],
          &data13[10 * map4idx],
          &data14[12 * map4idx]);
      }

    }
    OP_kernels[66].transfer  += Plan->transfer;
    OP_kernels[66].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[66].time     += wall_t2 - wall_t1;
}
