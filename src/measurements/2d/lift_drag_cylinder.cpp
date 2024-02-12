#include "measurements/2d/lift_drag_cylinder.h"

#include "op_seq.h"

#include "dg_dat_pool.h"
#include "dg_abort.h"

extern DGDatPool *dg_dat_pool;

LiftDragCylinder2D::LiftDragCylinder2D(SimulationDriver *d, const DG_FP refMu, const DG_FP x0, 
                                       const DG_FP y0, const DG_FP x1, const DG_FP y1, 
                                       const int sample_iter) : Measurement2D(d, sample_iter) {
  mu = refMu;
  box[0] = x0;
  box[1] = y0;
  box[2] = x1;
  box[3] = y1;
  if(dynamic_cast<INSSolverBase2D*>(d) == nullptr) {
    dg_abort("L2ErrorVortex2D measurement can only be used with 2D single or multiphase solver\n");
  }
  ins = dynamic_cast<INSSolverBase2D*>(d);
}

void LiftDragCylinder2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = ins->get_mesh();
  op_dat u = ins->get_vel_x();
  op_dat v = ins->get_vel_y();
  op_dat pr = ins->get_pr();

  DGTempDat ux = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat uy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat vx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat vy = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->grad(u, ux.dat, uy.dat);
  mesh->grad(v, vx.dat, vy.dat);

  DG_FP lift_coeff = 0.0;
  DG_FP drag_coeff = 0.0;
  if(mesh->bface2cells) {
    op_par_loop(measure_lift_drag, "measure_lift_drag", mesh->bfaces,
                op_arg_gbl(box, 4, DG_FP_STR, OP_READ),
                op_arg_gbl(&lift_coeff, 1, DG_FP_STR, OP_INC),
                op_arg_gbl(&drag_coeff, 1, DG_FP_STR, OP_INC),
                op_arg_gbl(&mu, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ,  -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(pr, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(ux.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(uy.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vx.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vy.dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ));
  }

  dg_dat_pool->releaseTempDatCells(ux);
  dg_dat_pool->releaseTempDatCells(uy);
  dg_dat_pool->releaseTempDatCells(vx);
  dg_dat_pool->releaseTempDatCells(vy);

  LiftDragHistory tmp;
  tmp.time = ins->get_time();
  tmp.lift = lift_coeff;
  tmp.drag = drag_coeff;
  history.push_back(tmp);
}

std::string LiftDragCylinder2D::get_filename() {
  return "lift_drag_coefficients";
}

std::string LiftDragCylinder2D::get_csv_header() {
  return "time,lift_coefficient,drag_coefficient";
}

std::string LiftDragCylinder2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].lift) + ",";
    result = result + double_to_text(history[io_count].drag);
    io_count++;
    return result;
  } else {
    return "";
  }
}