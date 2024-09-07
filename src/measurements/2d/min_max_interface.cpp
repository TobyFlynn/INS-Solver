#include "measurements/2d/min_max_interface.h"

#include "op_seq.h"

#include "dg_dat_pool.h"
#include "dg_abort.h"
#include "op2_utils.h"

extern DGDatPool *dg_dat_pool;

MinMaxInterface2D::MinMaxInterface2D(SimulationDriver *d, const int sample_iter) : Measurement2D(d, sample_iter) {
  if(dynamic_cast<MPINSSolver2D*>(d) == nullptr) {
    dg_abort("MinMaxInterface2D measurement can only be used with 2D multiphase solver\n");
  }

  mpins = dynamic_cast<MPINSSolver2D*>(d);
}

void MinMaxInterface2D::measure() {
  if(!sample_this_iter())
    return;

  DGMesh2D *mesh = mpins->get_mesh();
  op_dat ls = mpins->get_ls();
  LevelSetSolver2D *ls_solver = mpins->get_ls_solver();

  DG_FP min_x = 1e6;
  DG_FP min_y = 1e6;
  DG_FP max_x = -1e6;
  DG_FP max_y = -1e6;

  const DG_FP *s_ptr = getOP2PtrHostHE(ls, OP_READ);
  std::set<int> cellInds;
  for(int i = 0; i < mesh->cells->size; i++) {
    bool reinit = false;
    bool pos = s_ptr[i * DG_NP] >= 0.0;
    for(int j = 1; j < DG_NP; j++) {
      if((s_ptr[i * DG_NP + j] >= 0.0) != pos) {
        reinit = true;
      }
    }
    if(reinit) {
      cellInds.insert(i);
    }
  }

  const DG_FP *_x_ptr = getOP2PtrHostHE(mesh->x, OP_READ);
  const DG_FP *_y_ptr = getOP2PtrHostHE(mesh->y, OP_READ);

  std::map<int,int> _cell2polyMap;
  std::vector<PolyApprox> _polys;
  std::map<int,std::set<int>> stencils = PolyApprox::get_stencils(cellInds, mesh->face2cells, _x_ptr, _y_ptr);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    std::set<int> stencil = stencils.at(*it);
    PolyApprox p(*it, stencil, _x_ptr, _y_ptr, s_ptr, ls_solver->h);
    _polys.push_back(p);
    _cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostHE(mesh->x, OP_READ, _x_ptr);
  releaseOP2PtrHostHE(mesh->y, OP_READ, _y_ptr);
  releaseOP2PtrHostHE(ls, OP_READ, s_ptr);

  DGTempDat tmp_sampleX = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleY = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  ls_solver->sampleInterface(tmp_sampleX.dat, tmp_sampleY.dat, _polys, _cell2polyMap, cellInds);

  const DG_FP *sample_pts_x = getOP2PtrHost(tmp_sampleX.dat, OP_READ);
  const DG_FP *sample_pts_y = getOP2PtrHost(tmp_sampleY.dat, OP_READ);

  // Get min max from sample points
  for(int i = 0; i < mesh->cells->size * LS_SAMPLE_NP; i++) {
    if(!isnan(sample_pts_x[i]) && !isnan(sample_pts_y[i])) {
      min_x = sample_pts_x[i] < min_x ? sample_pts_x[i] : min_x;
      max_x = sample_pts_x[i] > max_x ? sample_pts_x[i] : max_x;
      min_y = sample_pts_y[i] < min_y ? sample_pts_y[i] : min_y;
      max_y = sample_pts_y[i] > max_y ? sample_pts_y[i] : max_y;
    }
  }

  releaseOP2PtrHost(tmp_sampleX.dat, OP_READ, sample_pts_x);
  releaseOP2PtrHost(tmp_sampleY.dat, OP_READ, sample_pts_y);
  dg_dat_pool->releaseTempDatCells(tmp_sampleX);
  dg_dat_pool->releaseTempDatCells(tmp_sampleY);

/*  
  op_par_loop(measure_min_max_interface, "measure_min_max_interface", mesh->cells,
              op_arg_gbl(&min_x, 1, DG_FP_STR, OP_MIN),
              op_arg_gbl(&min_y, 1, DG_FP_STR, OP_MIN),
              op_arg_gbl(&max_x, 1, DG_FP_STR, OP_MAX),
              op_arg_gbl(&max_y, 1, DG_FP_STR, OP_MAX),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(ls, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
*/

  DG_FP time = mpins->get_time();
  MinMaxHistory tmp;
  tmp.time = time;
  tmp.min_x = min_x;
  tmp.min_y = min_y;
  tmp.max_x = max_x;
  tmp.max_y = max_y;
  history.push_back(tmp);
}

std::string MinMaxInterface2D::get_filename() {
  return "min_max_interface";
}

std::string MinMaxInterface2D::get_csv_header() {
  return "time,min_x,min_y,max_x,max_y";
}

std::string MinMaxInterface2D::get_next_csv_line() {
  if(io_count < history.size()) {
    std::string result = double_to_text(history[io_count].time) + ",";
    result = result + double_to_text(history[io_count].min_x) + ",";
    result = result + double_to_text(history[io_count].min_y) + ",";
    result = result + double_to_text(history[io_count].max_x) + ",";
    result = result + double_to_text(history[io_count].max_y);
    io_count++;
    return result;
  } else {
    return "";
  }
}