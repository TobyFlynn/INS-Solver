//
// auto-generated by op2.py
//

//user function
#include "../kernels/update_Q.h"

// host stub function
void op_par_loop_update_Q(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(67);
  OP_kernels[67].name      = name;
  OP_kernels[67].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  update_Q");
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
        update_Q(
          (double*)arg0.data,
          &((double*)arg1.data)[15*n],
          &((double*)arg2.data)[15*n],
          &((double*)arg3.data)[15*n],
          &((double*)arg4.data)[15*n]);
      }
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[67].time     += wall_t2 - wall_t1;
  OP_kernels[67].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[67].transfer += (float)set->size * arg2.size;
  OP_kernels[67].transfer += (float)set->size * arg3.size;
  OP_kernels[67].transfer += (float)set->size * arg4.size;
}
