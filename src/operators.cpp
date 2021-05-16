#include "op_seq.h"
#include "ins_data.h"
#include "blas_calls.h"

#include "kernels/div.h"
#include "kernels/curl.h"
#include "kernels/grad.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res) {
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DR), 15, u, 0.0, data->div[0]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DS), 15, u, 0.0, data->div[1]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DR), 15, v, 0.0, data->div[2]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DS), 15, v, 0.0, data->div[3]);

  op_par_loop(div, "div", data->cells,
              op_arg_dat(data->div[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ry, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(res, -1, OP_ID, 15, "double", OP_WRITE));
}

void curl(INSData *data, op_dat u, op_dat v, op_dat res) {
  // Same matrix multiplications as div
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DR), 15, u, 0.0, data->div[0]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DS), 15, u, 0.0, data->div[1]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DR), 15, v, 0.0, data->div[2]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DS), 15, v, 0.0, data->div[3]);

  op_par_loop(curl, "curl", data->cells,
              op_arg_dat(data->div[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ry, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(res, -1, OP_ID, 15, "double", OP_WRITE));
}

void grad(INSData *data, op_dat u, op_dat ux, op_dat uy) {
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DR), 15, u, 0.0, data->div[0]);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DS), 15, u, 0.0, data->div[1]);

  op_par_loop(grad, "grad", data->cells,
              op_arg_dat(data->div[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ry, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(ux, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(uy, -1, OP_ID, 15, "double", OP_WRITE));
}
