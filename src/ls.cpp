#include "ls.h"

#include "op_seq.h"

#include "kernels/init_surface.h"

LS::LS(INSData *d) {
  data = d;

  s_data = (double *)calloc(15 * data->numCells, sizeof(double));

  rk_data[0] = (double *)calloc(15 * data->numCells, sizeof(double));
  rk_data[1] = (double *)calloc(15 * data->numCells, sizeof(double));
  rk_data[2] = (double *)calloc(15 * data->numCells, sizeof(double));
  rkQ_data   = (double *)calloc(15 * data->numCells, sizeof(double));

  s = op_decl_dat(data->cells, 15, "double", s_data, "s");

  rk[0] = op_decl_dat(data->cells, 15, "double", rk_data[0], "rk0");
  rk[1] = op_decl_dat(data->cells, 15, "double", rk_data[1], "rk1");
  rk[2] = op_decl_dat(data->cells, 15, "double", rk_data[2], "rk2");
  rkQ   = op_decl_dat(data->cells, 15, "double", rkQ_data, "rkQ");
}

LS::~LS() {
  free(s_data);

  free(rk_data[0]);
  free(rk_data[1]);
  free(rk_data[2]);
  free(rkQ_data);
}

void LS::init() {
  op_par_loop(init_surface, "init_surface", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, 15, "double", OP_WRITE));
}
