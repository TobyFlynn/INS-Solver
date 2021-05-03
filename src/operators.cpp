#include "op_seq.h"
#include "ins_data.h"
#include "blas_calls.h"

#include "kernels/div.h"
#include "kernels/curl.h"
#include "kernels/grad.h"
#include "kernels/cub_grad.h"
#include "kernels/cub_div.h"

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

void cub_grad(INSData *data, CubatureData *cubData, op_dat u, op_dat ux, op_dat uy, bool weak) {
  if(weak) {
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, u, 0.0, cubData->op_temps[0]);

    op_par_loop(cub_grad_w, "cub_grad_w", data->cells,
                op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
                op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_WRITE),
                op_arg_dat(cubData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
                op_arg_dat(cubData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cubData->op_temps[0], 0.0, ux);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cubData->op_temps[1], 1.0, ux);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cubData->op_temps[2], 0.0, uy);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cubData->op_temps[3], 1.0, uy);
  } else {
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, u, 0.0, cubData->op_temps[0]);
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, u, 0.0, cubData->op_temps[1]);

    op_par_loop(cub_grad, "cub_grad", data->cells,
                op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
                op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_RW));

    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, cubData->op_temps[0], 0.0, ux);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, cubData->op_temps[1], 0.0, uy);
  }
}

void cub_div(INSData *data, CubatureData *cubData, op_dat u, op_dat v, op_dat res, bool weak) {
  if(weak) {
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, u, 0.0, cubData->op_temps[0]);
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, v, 0.0, cubData->op_temps[1]);

    op_par_loop(cub_div_w, "cub_div_w", data->cells,
                op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
                op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_RW),
                op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
                op_arg_dat(cubData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cubData->op_temps[0], 0.0, res);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cubData->op_temps[1], 1.0, res);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cubData->op_temps[2], 1.0, res);
    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cubData->op_temps[3], 1.0, res);
  } else {
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, u, 0.0, cubData->op_temps[0]);
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, u, 0.0, cubData->op_temps[1]);
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, v, 0.0, cubData->op_temps[2]);
    op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, v, 0.0, cubData->op_temps[3]);

    op_par_loop(cub_div, "cub_div", data->cells,
                op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->J, -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
                op_arg_dat(cubData->op_temps[1], -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[2], -1, OP_ID, 46, "double", OP_READ),
                op_arg_dat(cubData->op_temps[3], -1, OP_ID, 46, "double", OP_READ));

    op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, cubData->op_temps[0], 0.0, res);
  }
}
