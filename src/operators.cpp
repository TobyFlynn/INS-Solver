#include "op_seq.h"
#include "ins_data.h"
#include "blas_calls.h"

#include "kernels/div.h"
#include "kernels/curl.h"
#include "kernels/grad.h"
#include "kernels/cub_grad.h"
#include "kernels/cub_div.h"

void div(INSData *data, op_dat u, op_dat v, op_dat res) {
  div_blas(data, u, v);

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
  div_blas(data, u, v);

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
  grad_blas(data, u);

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

void cub_grad(INSData *data, CubatureData *cubData, op_dat u, op_dat ux, op_dat uy) {
  cub_grad_blas(data, cubData, u);

  op_par_loop(cub_grad, "cub_grad", data->cells,
              op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cubData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cubData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

  cub_grad_blas2(data, cubData, ux, uy);
}

void cub_div(INSData *data, CubatureData *cubData, op_dat u, op_dat v, op_dat res) {
  cub_div_blas(data, cubData, u, v);

  op_par_loop(cub_div, "cub_div", data->cells,
              op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cubData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cubData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

  cub_div_blas2(data, cubData, res);
}
