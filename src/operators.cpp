#include "operators.h"

#include "op_seq.h"
#include "ins_data.h"
#include "blas_calls.h"

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

void cub_grad(INSData *data, CubatureData *cData, op_dat u, op_dat ux, op_dat uy) {
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DR), 15, u, 0.0, cData->op_temps[0]);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_DS), 15, u, 0.0, cData->op_temps[1]);

  op_par_loop(cub_grad, "cub_grad", data->cells,
              op_arg_dat(cData->rx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->ry, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sy, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cData->op_temps[1], -1, OP_ID, 46, "double", OP_RW));

  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, cData->op_temps[0], 0.0, ux);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_V), 15, cData->op_temps[1], 0.0, uy);

  inv_mass(data, ux);
  inv_mass(data, uy);
}

void cub_grad_weak(INSData *data, CubatureData *cData, op_dat u, op_dat ux, op_dat uy) {
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, u, 0.0, cData->op_temps[0]);

  op_par_loop(cub_grad_weak, "cub_grad_weak", data->cells,
              op_arg_dat(cData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cData->rx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->ry, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sy, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->op_temps[1], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cData->op_temps[0], 0.0, ux);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cData->op_temps[1], 1.0, ux);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cData->op_temps[2], 0.0, uy);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cData->op_temps[3], 1.0, uy);

  // inv_mass(data, ux);
  // inv_mass(data, uy);
}

void cub_div_weak(INSData *data, CubatureData *cData, op_dat u, op_dat v, op_dat res) {
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, u, 0.0, cData->op_temps[0]);
  op2_gemv(true, 46, 15, 1.0, constants->get_ptr(Constants::CUB_V), 15, v, 0.0, cData->op_temps[1]);

  op_par_loop(cub_div_weak, "cub_div_weak", data->cells,
              op_arg_dat(cData->op_temps[0], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cData->op_temps[1], -1, OP_ID, 46, "double", OP_RW),
              op_arg_dat(cData->rx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sx, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->ry, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->sy, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->J, -1, OP_ID, 46, "double", OP_READ),
              op_arg_dat(cData->op_temps[2], -1, OP_ID, 46, "double", OP_WRITE),
              op_arg_dat(cData->op_temps[3], -1, OP_ID, 46, "double", OP_WRITE));

  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cData->op_temps[0], 0.0, res);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cData->op_temps[1], 1.0, res);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DR), 15, cData->op_temps[2], 1.0, res);
  op2_gemv(false, 15, 46, 1.0, constants->get_ptr(Constants::CUB_DS), 15, cData->op_temps[3], 1.0, res);
}

void inv_mass(INSData *data, op_dat u) {
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::INV_MASS), 15, u, 0.0, data->div[0]);

  op_par_loop(inv_J, "inv_J", data->cells,
              op_arg_dat(data->J, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->div[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE));
}
