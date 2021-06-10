//
// auto-generated by op2.py
//

//user function
#include "../kernels/viscosity_reset_bc.h"

// host stub function
void op_par_loop_viscosity_reset_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(52);
  OP_kernels[52].name      = name;
  OP_kernels[52].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_reset_bc");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  if (set_size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        viscosity_reset_bc(
          &((double*)arg0.data)[21*n],
          &((double*)arg1.data)[21*n]);
      }
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[52].time     += wall_t2 - wall_t1;
  OP_kernels[52].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[52].transfer += (float)set->size * arg1.size * 2.0f;
}
