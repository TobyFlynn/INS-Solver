//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void save_values_openacc( const double *v_vals, double *c_vals) {
  #if DG_ORDER == 4
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[5]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[5] + v_vals[6]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[6]) / 3.0;
  c_vals[3]  = (v_vals[2] + v_vals[6] + v_vals[7]) / 3.0;
  c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[7]) / 3.0;
  c_vals[5]  = (v_vals[3] + v_vals[7] + v_vals[8]) / 3.0;
  c_vals[6]  = (v_vals[3] + v_vals[4] + v_vals[8]) / 3.0;
  c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[9]) / 3.0;
  c_vals[8]  = (v_vals[6] + v_vals[9] + v_vals[10]) / 3.0;
  c_vals[9]  = (v_vals[6] + v_vals[7] + v_vals[10]) / 3.0;
  c_vals[10] = (v_vals[7] + v_vals[10] + v_vals[11]) / 3.0;
  c_vals[11] = (v_vals[7] + v_vals[8] + v_vals[11]) / 3.0;
  c_vals[12] = (v_vals[9] + v_vals[10] + v_vals[12]) / 3.0;
  c_vals[13] = (v_vals[10] + v_vals[12] + v_vals[13]) / 3.0;
  c_vals[14] = (v_vals[10] + v_vals[11] + v_vals[13]) / 3.0;
  c_vals[15] = (v_vals[12] + v_vals[13] + v_vals[14]) / 3.0;
  #elif DG_ORDER == 3
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[4]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[4] + v_vals[5]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[5]) / 3.0;
  c_vals[3]  = (v_vals[2] + v_vals[5] + v_vals[6]) / 3.0;
  c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[6]) / 3.0;
  c_vals[5]  = (v_vals[4] + v_vals[5] + v_vals[7]) / 3.0;
  c_vals[6]  = (v_vals[5] + v_vals[7] + v_vals[8]) / 3.0;
  c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[8]) / 3.0;
  c_vals[8]  = (v_vals[7] + v_vals[8] + v_vals[9]) / 3.0;
  #elif DG_ORDER == 2
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[3]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[3] + v_vals[4]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[4]) / 3.0;
  c_vals[3]  = (v_vals[3] + v_vals[4] + v_vals[5]) / 3.0;
  #else
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[2]) / 3.0;
  #endif
}

// host stub function
void op_par_loop_save_values(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(48);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[48].name      = name;
  OP_kernels[48].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  save_values");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1)
    for ( int n=0; n<set->size; n++ ){
      save_values_openacc(
        &data0[10*n],
        &data1[9*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[48].time     += wall_t2 - wall_t1;
  OP_kernels[48].transfer += (float)set->size * arg0.size;
  OP_kernels[48].transfer += (float)set->size * arg1.size * 2.0f;
}