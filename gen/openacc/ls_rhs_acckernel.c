//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_rhs_openacc( const double *sign, const double *dpldx, const double *dprdx,
                   const double *dpldy, const double *dprdy, double *rk) {
  for(int i = 0; i < 10; i++) {
    if(sign[i] > 0.0) {
      double plmx2 = fmin(dpldx[i], 0.0);
      plmx2 = plmx2 * plmx2;
      double plmy2 = fmin(dpldy[i], 0.0);
      plmy2 = plmy2 * plmy2;
      double prpx2 = fmax(dprdx[i], 0.0);
      prpx2 = prpx2 * prpx2;
      double prpy2 = fmax(dprdy[i], 0.0);
      prpy2 = prpy2 * prpy2;

      double H = sign[i] * (sqrt(fmax(plmx2, prpx2) + fmax(plmy2, prpy2)) - 1.0);
      rk[i] = -H;
    } else {
      double plpx2 = fmax(dpldx[i], 0.0);
      plpx2 = plpx2 * plpx2;
      double plpy2 = fmax(dpldy[i], 0.0);
      plpy2 = plpy2 * plpy2;
      double prmx2 = fmin(dprdx[i], 0.0);
      prmx2 = prmx2 * prmx2;
      double prmy2 = fmin(dprdy[i], 0.0);
      prmy2 = prmy2 * prmy2;

      double H = sign[i] * (sqrt(fmax(plpx2, prmx2) + fmax(plpy2, prmy2)) - 1.0);
      rk[i] = -H;
    }
  }
}

// host stub function
void op_par_loop_ls_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(58);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[58].name      = name;
  OP_kernels[58].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      ls_rhs_openacc(
        &data0[10*n],
        &data1[10*n],
        &data2[10*n],
        &data3[10*n],
        &data4[10*n],
        &data5[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[58].time     += wall_t2 - wall_t1;
  OP_kernels[58].transfer += (float)set->size * arg0.size;
  OP_kernels[58].transfer += (float)set->size * arg1.size;
  OP_kernels[58].transfer += (float)set->size * arg2.size;
  OP_kernels[58].transfer += (float)set->size * arg3.size;
  OP_kernels[58].transfer += (float)set->size * arg4.size;
  OP_kernels[58].transfer += (float)set->size * arg5.size * 2.0f;
}
