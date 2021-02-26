//
// auto-generated by op2.py
//

#include  "op_lib_cpp.h"

//
// op_par_loop declarations
//
#ifdef OPENACC
#ifdef __cplusplus
extern "C" {
#endif
#endif

void op_par_loop_div(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_curl(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_grad(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );
#ifdef OPENACC
#ifdef __cplusplus
}
#endif
#endif

#include "ins_data.h"
#include "blas_calls.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res, bool weak) {
  div_blas(data, u, v, weak);

  op_par_loop_div("div",data->cells,
              op_arg_dat(data->div[0],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[1],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[2],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[3],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->rx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->ry,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sy,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(res,-1,OP_ID,15,"double",OP_WRITE));
}

void curl(INSData *data, op_dat u, op_dat v, op_dat res) {
  // Same matrix multiplications as div
  // Rename this later
  div_blas(data, u, v);

  op_par_loop_curl("curl",data->cells,
              op_arg_dat(data->div[0],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[1],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[2],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[3],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->rx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->ry,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sy,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(res,-1,OP_ID,15,"double",OP_WRITE));
}

void grad(INSData *data, op_dat u, op_dat ux, op_dat uy, bool weak) {
  grad_blas(data, u, weak);

  op_par_loop_grad("grad",data->cells,
              op_arg_dat(data->div[0],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->div[1],-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->rx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sx,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->ry,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(data->sy,-1,OP_ID,15,"double",OP_READ),
              op_arg_dat(ux,-1,OP_ID,15,"double",OP_WRITE),
              op_arg_dat(uy,-1,OP_ID,15,"double",OP_WRITE));
}
