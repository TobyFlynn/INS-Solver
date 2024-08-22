#include "solvers/2d/ls_solver.h"

#include "op_seq.h"

#include <limits>
#include <cmath>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "timing.h"
#include "config.h"
#include "op2_utils.h"
#ifdef INS_MPI
#include "ls_utils/2d/kd_tree_mpi.h"
#endif

extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;
extern Timing *timer;
extern Config *config;

/**************************************************************************
 * LS Advection Solver class that extends the base Advection Solver class *
 **************************************************************************/
LevelSetAdvectionSolver2D::LevelSetAdvectionSolver2D(DGMesh2D *m) : AdvectionSolver2D(m) {}

void LevelSetAdvectionSolver2D::set_bc_types(op_dat bc) {
  bc_types = bc;
}

void LevelSetAdvectionSolver2D::bc_kernel(op_dat val, op_dat u, op_dat v, op_dat out) {
  op_par_loop(ls_advec_2d_bc, "ls_advec_2d_bc", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
}

void LevelSetAdvectionSolver2D::bc_kernel_oi(op_dat val, op_dat u, op_dat v, op_dat flux) {
  op_par_loop(ls_advec_2d_oi_bc, "ls_advec_2d_oi_bc", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));
}

/************************
 * Main LS Solver class *
 ************************/
LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m) {
  mesh = m;
  resuming = false;

  s = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_s");
  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdy");
  kink = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_kink");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_y");

  advecSolver = new LevelSetAdvectionSolver2D(mesh);
}

LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m, const std::string &filename) {
  mesh = m;
  resuming = true;

  s = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ls_s");

  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdy");
  kink = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_kink");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_y");

  advecSolver = new LevelSetAdvectionSolver2D(mesh);
}

LevelSetSolver2D::~LevelSetSolver2D() {
  delete advecSolver;
}

void LevelSetSolver2D::init() {
  if(!resuming) {
    op_par_loop(init_surface_2d, "init_surface_2d", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(s,       -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  // h = std::numeric_limits<DG_FP>::max();
  h = 0.0;
  op_par_loop(calc_h_ls, "calc_h_ls", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  
  op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
              op_arg_dat(kink, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op_printf("LS h: %g\n", h);
  // alpha = 2.0 * h / DG_ORDER;
  // order_width = 2.0 * h;
  // epsilon = h / DG_ORDER;
  alpha = 12.0 * h;
  // order_width = 12.0 * h;
  epsilon = h;
  // reinit_width = 20.0 * h;
  ls_cap = 50.0;
  config->getDouble("level-set-options", "ls_cap", ls_cap);
  reinit_width = ls_cap;
  // reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  // numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  reinit_counter = 0;
  reinit_frequency = 14;
  config->getInt("level-set-options", "reinit_freq", reinit_frequency);
  reinitialise = reinit_frequency > 0;

  op_printf("Alpha: %g\t\tReinit Width: %g\n", alpha, reinit_width);

  #ifdef INS_MPI
  kdtree = new KDTreeMPI(mesh, 2.0 * reinit_width);
  #else
  kdtree = new KDTree(mesh);
  #endif

  int tmp_reinit_ic = 1;
  config->getInt("level-set-options", "reinit_ic", tmp_reinit_ic);
  
  int tmp_kink = 1;
  config->getInt("level-set-options", "kink_detection", tmp_kink);
  kink_detection = tmp_kink == 1;
  int tmp_kink_element = 0;
  config->getInt("level-set-options", "kink_avoid_whole_element", tmp_kink_element);
  kink_avoid_whole_element = tmp_kink_element == 1;

  kink_max_distance_between_points = 2.0;
  config->getDouble("level-set-options", "kink_max_distance_between_points", kink_max_distance_between_points);
  kink_max_distance_between_points *= h;
  kink_max_neighbours = 100;
  config->getInt("level-set-options", "kink_max_neighbours", kink_max_neighbours);
  kink_sqr_tol = 0.5;
  config->getDouble("level-set-options", "kink_tol", kink_sqr_tol);
  kink_sqr_tol *= kink_sqr_tol;
  if(kink_detection)
    create_point_map_for_kink_detection();

  if(tmp_reinit_ic != 0)
    reinitLS();
}

void LevelSetSolver2D::set_bc_types(op_dat bc) {
  advecSolver->set_bc_types(bc);
}

void LevelSetSolver2D::setVelField(op_dat u1, op_dat v1) {
  u = u1;
  v = v1;
}

void LevelSetSolver2D::step(const DG_FP dt, const int num_steps) {
  timer->startTimer("LevelSetSolver2D - step");
  advecSolver->set_dt(dt);
  for(int i = 0; i < num_steps; i++)
    advecSolver->step(s, u, v);

  reinit_counter++;
  if(reinitialise && reinit_counter >= reinit_frequency) {
    reinitLS();
    reinit_counter = 0;
  }
  timer->endTimer("LevelSetSolver2D - step");
}

bool LevelSetSolver2D::reinitNeeded() {
  DG_FP res = 0.0;
  int count = 0;
  mesh->grad(s, dsdx, dsdy);
  op_par_loop(ls_reinit_check, "ls_reinit_check", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dsdx,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dsdy,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&res,   1, DG_FP_STR, OP_INC),
              op_arg_gbl(&count, 1, "int", OP_INC));

  res = res / (DG_FP)count;
  // std::cout << "LS residual: " << res << " " << abs(1.0 - res) << std::endl;
  return abs(1.0 - res) > 0.01;
}

void LevelSetSolver2D::getRhoMu(op_dat rho, op_dat mu) {
  timer->startTimer("LevelSetSolver2D - getRhoMu");
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver2D - getRhoMu");
}

void LevelSetSolver2D::getRhoVolOI(op_dat rho) {
  timer->startTimer("LevelSetSolver2D - getRhoVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, s, 0.0, rho);
  op_par_loop(ls_step_rho_vol_oi, "ls_step_rho_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getRhoVolOI");
}

void LevelSetSolver2D::getMuVolOI(op_dat mu) {
  timer->startTimer("LevelSetSolver2D - getMuVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, s, 0.0, mu);
  op_par_loop(ls_step_mu_vol_oi, "ls_step_mu_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(mu, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getMuVolOI");
}

void LevelSetSolver2D::getRhoSurfOI(op_dat rho) {
  timer->startTimer("LevelSetSolver2D - getRhoSurfOI");
  DGTempDat rho_tmp = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  op_par_loop(ls_step_surf_oi_0, "ls_step_surf_oi_0", mesh->cells,
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho_tmp.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, rho_tmp.dat, 0.0, rho);
  dg_dat_pool->releaseTempDatCells(rho_tmp);
  op_par_loop(ls_step_surf_oi_1, "ls_step_surf_oi_1", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getRhoSurfOI");
}

void LevelSetSolver2D::getNormalsCurvature(op_dat nx, op_dat ny, op_dat curv) {
  timer->startTimer("LevelSetSolver2D - getNormalsCurvature");
  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  mesh->grad(s, nx, ny);
  mesh->div(nx, ny, curv);
  timer->endTimer("LevelSetSolver2D - getNormalsCurvature");
}

bool newtoncp_gepp(arma::mat &A, arma::vec &b) {
  for (int i = 0; i < 3; ++i) {
    int j = i;
    for (int k = i + 1; k < 3; ++k)
      if (std::abs(A(k,i)) > std::abs(A(j,i)))
        j = k;
    if (j != i) {
      for (int k = 0; k < 3; ++k)
        std::swap(A(i,k), A(j,k));
      std::swap(b(i), b(j));
    }

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<DG_FP>::epsilon())
      return false;

    DG_FP fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 3; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    DG_FP sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(DG_FP &closest_pt_x, DG_FP &closest_pt_y, const DG_FP node_x,
                   const DG_FP node_y, PolyApprox &p, const DG_FP h) {
  DG_FP lambda = 0.0;
  bool converged = false;
  DG_FP pt_x = closest_pt_x;
  DG_FP pt_y = closest_pt_y;
  DG_FP init_x = closest_pt_x;
  DG_FP init_y = closest_pt_y;

  for(int step = 0; step < 50; step++) {
    DG_FP pt_x_old = pt_x;
    DG_FP pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    DG_FP surface = p.val_at(pt_x, pt_y);
    DG_FP surface_dx, surface_dy;
    p.grad_at(pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    DG_FP hessian[3];
    p.hessian_at(pt_x, pt_y, hessian[0], hessian[1], hessian[2]);

    // Check if |nabla(surface)| = 0, if so then return
    DG_FP gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;

    // Gradient of functional
    arma::vec gradf(3);
    gradf(0) = pt_x - node_x + lambda * surface_dx;
    gradf(1) = pt_y - node_y + lambda * surface_dy;
    gradf(2) = surface;

    // Calculate Hessian of functional
    arma::mat hessianf(3, 3);
    hessianf(0, 0) = 1.0 + lambda * hessian[0];
    hessianf(0, 1) = lambda * hessian[1]; hessianf(1, 0) = hessianf(0, 1);
    hessianf(0, 2) = surface_dx; hessianf(2, 0) = hessianf(0, 2);

    hessianf(1, 1) = 1.0 + lambda * hessian[2];
    hessianf(1, 2) = surface_dy; hessianf(2, 1) = hessianf(1, 2);

    hessianf(2, 2) = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      DG_FP delta1_x = (surface / gradsqrnorm) * surface_dx;
      DG_FP delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;
      DG_FP delta2_x = pt_x - node_x + lambda * surface_dx;
      DG_FP delta2_y = pt_y - node_y + lambda * surface_dy;
      DG_FP msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      DG_FP msqr = ans(0) * ans(0) + ans(1) * ans(1);
      if(msqr > h * 0.5 * h * 0.5)
        ans = ans * 0.5 * h / sqrt(msqr);

      // Update guess
      pt_x -= ans(0);
      pt_y -= ans(1);
      lambda -= ans(2);
    }

    // Gone outside the element, return
    if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > 1.5 * 1.5 * h * h) {
      pt_x = pt_x_old;
      pt_y = pt_y_old;
      break;
    }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-24) {
      converged = true;
      break;
    }
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;

  return converged;
}

void newton_method(const int numPts, DG_FP *closest_x, DG_FP *closest_y,
                   const DG_FP *x, const DG_FP *y, int *poly_ind,
                   std::vector<PolyApprox> &polys, DG_FP *s, DG_FP *kink, const DG_FP h,
                   const DG_FP reinit_width, const DG_FP ls_cap, bool kink_avoid_whole_element) {
  int numNonConv = 0;
  int numReinit = 0;

  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    DG_FP dist2 = (closest_x[start_ind] - x[start_ind]) * (closest_x[start_ind] - x[start_ind])
                + (closest_y[start_ind] - y[start_ind]) * (closest_y[start_ind] - y[start_ind]);

    if(dist2 < (reinit_width + 2.0 * h) * (reinit_width + 2.0 * h)) {
      DG_FP dist1 = (closest_x[i] - x[i]) * (closest_x[i] - x[i])
                  + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
      
      if(kink_avoid_whole_element) {
        bool kink_cell = false;
        for(int n = 0; n < DG_NP; n++) {
          kink_cell = kink_cell || (kink[start_ind + n] != 0.0);
        }
        if(kink_cell) continue;
      } else {
        if(kink[i] != 0.0 && dist1 < 9.0 * h * h)
          continue;
        else
          kink[i] = 0.0;
      }

      DG_FP poly_offset_x, poly_offset_y;
      polys[poly_ind[i]].get_offsets(poly_offset_x, poly_offset_y);
      DG_FP _closest_x = closest_x[i] - poly_offset_x;
      DG_FP _closest_y = closest_y[i] - poly_offset_y;
      DG_FP _node_x = x[i] - poly_offset_x;
      DG_FP _node_y = y[i] - poly_offset_y;
      bool converged = newton_kernel(_closest_x, _closest_y, _node_x, _node_y, polys[poly_ind[i]], h);
      if(converged) {
        // DG_FP dsdx, dsdy;
        // polys[poly_ind[i]].grad_at(_closest_x, _closest_y, dsdx, dsdy);
        // DG_FP dot = (_node_x - _closest_x) * dsdx + (_node_y - _closest_y) * dsdy;
        // bool negative = dot < 0.0;
        bool negative = s[i] < 0.0;
        s[i] = (_closest_x - _node_x) * (_closest_x - _node_x) + (_closest_y - _node_y) * (_closest_y - _node_y);
        s[i] = sqrt(s[i]);
        if(negative) s[i] *= -1.0;
      } else {
        bool negative = s[i] < 0.0;
        // s[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
        s[i] = dist1;
        s[i] = sqrt(s[i]);
        if(negative) s[i] *= -1.0;
        #pragma omp atomic
        numNonConv++;
      }
      #pragma omp atomic
      numReinit++;
    } else {
      bool negative = s[i] < 0.0;
      s[i] = negative ? -ls_cap : ls_cap;
    }
  }

  double percent_non_converge = numReinit == 0 ? 0.0 : (double)numNonConv / (double)numReinit;
  if(percent_non_converge > 0.1)
    std::cout << percent_non_converge * 100.0 << "% reinitialisation points did not converge" << std::endl;
}

struct LSKinkPoint {
  DGUtils::Vec<2> coord;
  DG_FP s, dsdx, dsdy;
  int count;
  std::vector<int> indices;
};

void LevelSetSolver2D::create_point_map_for_kink_detection() {
  timer->startTimer("LevelSetSolver2D - point map creation");
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);

  std::map<DGUtils::Vec<2>, LSKinkPoint> pointMap;
  for(int i = 0; i < mesh->cells->size * DG_NP; i++) {
    DGUtils::Vec<2> coord(x_ptr[i], y_ptr[i]);
    LSKinkPoint point;
    auto res = pointMap.insert({coord, point});
    if(res.second) {
      // Point was inserted
      res.first->second.coord = coord;
      res.first->second.count = 1;
      res.first->second.indices.push_back(i);
    } else {
      // Point already exists
      res.first->second.count++;
      res.first->second.indices.push_back(i);
    }
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);

  // Create a map consisting of buckets that contain a vector of points.
  // Each bucket is a square of length max_distance_between_points.
  // Key to the bucket is the coordinates of its bottom left vertex.
  const DG_FP max_distance_between_points_sqr = kink_max_distance_between_points * kink_max_distance_between_points;
  std::map<DGUtils::Vec<2>, std::vector<LSKinkPoint>> bucket_map;
  for(auto const &p : pointMap) {
    LSKinkPoint point = p.second;

    DGUtils::Vec<2> bucket_coord = point.coord;
    bucket_coord[0] = static_cast<int>(bucket_coord[0] / kink_max_distance_between_points) * kink_max_distance_between_points;
    bucket_coord[1] = static_cast<int>(bucket_coord[1] / kink_max_distance_between_points) * kink_max_distance_between_points;

    auto res = bucket_map.insert({bucket_coord, std::vector<LSKinkPoint>()});
    res.first->second.push_back(point);
  }

  DGUtils::Vec<2> bucket_offsets[8] = {
    DGUtils::Vec<2>(-kink_max_distance_between_points, 0.0),
    DGUtils::Vec<2>(-kink_max_distance_between_points, kink_max_distance_between_points),
    DGUtils::Vec<2>(0.0, kink_max_distance_between_points),
    DGUtils::Vec<2>(kink_max_distance_between_points, kink_max_distance_between_points),
    DGUtils::Vec<2>(kink_max_distance_between_points, 0.0),
    DGUtils::Vec<2>(kink_max_distance_between_points, -kink_max_distance_between_points),
    DGUtils::Vec<2>(0.0, -kink_max_distance_between_points),
    DGUtils::Vec<2>(-kink_max_distance_between_points, -kink_max_distance_between_points)
  };

  for(auto const &bucket : bucket_map) {
    // Get points in this bucket
    std::vector<LSKinkPoint> points(bucket.second);
    // Get points in surrounding buckets
    for(int i = 0; i < 8; i++) {
      DGUtils::Vec<2> bucket_coord = bucket.first + bucket_offsets[i];
      auto const neighbouring_bucket = bucket_map.find(bucket_coord);
      if(neighbouring_bucket != bucket_map.end()) {
        points.insert(points.end(), neighbouring_bucket->second.begin(), neighbouring_bucket->second.end());
      }
    }

    // Iterate over all points
    for(int i = 0; i < points.size(); i++) {
      std::vector<DGUtils::Vec<2>> pts;
      for(int j = 0; j < points.size(); j++) {
        if(j == i || (points[i].coord - points[j].coord).sqr_magnitude() > max_distance_between_points_sqr)
          continue;
        
        if(pts.size() < kink_max_neighbours) {
          pts.push_back(points[j].coord);
        } else {
          double sqr_dist_j = (points[i].coord - points[j].coord).sqr_magnitude();
          double max_distance = 0.0;
          int ind = -1;
          for(int k = 0; k < pts.size(); k++) {
            double sqr_dist_k = (points[i].coord - points[k].coord).sqr_magnitude();
            if(sqr_dist_j < sqr_dist_k) {
              if(sqr_dist_k > max_distance) {
                max_distance = sqr_dist_k;
                ind = k;
              }
            }
          }
          if(ind >= 0)
            pts[ind] = points[j].coord;
        }
      }
      point_map_for_kink_detection.insert({points[i].coord, pts});
    }
  }
  timer->endTimer("LevelSetSolver2D - point map creation");
}

void LevelSetSolver2D::detect_kinks(std::set<int> &stencilInds) {
  op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
              op_arg_dat(kink, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat dsdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dsdy = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(s, dsdx.dat, dsdy.dat);

  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  const DG_FP *dsdx_ptr = getOP2PtrHost(dsdx.dat, OP_READ);
  const DG_FP *dsdy_ptr = getOP2PtrHost(dsdy.dat, OP_READ);

  std::map<DGUtils::Vec<2>, LSKinkPoint> pointMap;

  for(int i = 0; i < mesh->cells->size * DG_NP; i++) {
    // if(fabs(s_ptr[i]) < 4.0 * h) {
    if(stencilInds.count(i / DG_NP) != 0) {
      DGUtils::Vec<2> coord(x_ptr[i], y_ptr[i]);
      LSKinkPoint point;
      auto res = pointMap.insert({coord, point});
      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.s = s_ptr[i];
        res.first->second.dsdx = dsdx_ptr[i];
        res.first->second.dsdy = dsdy_ptr[i];
        res.first->second.count = 1;
        res.first->second.indices.push_back(i);
      } else {
        // Point already exists
        res.first->second.s += s_ptr[i];
        res.first->second.dsdx += dsdx_ptr[i];
        res.first->second.dsdy += dsdy_ptr[i];
        res.first->second.count++;
        res.first->second.indices.push_back(i);
      }
    }
  }

  releaseOP2PtrHost(dsdy.dat, OP_READ, dsdy_ptr);
  releaseOP2PtrHost(dsdx.dat, OP_READ, dsdx_ptr);
  releaseOP2PtrHost(s, OP_READ, s_ptr);
  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);

  dg_dat_pool->releaseTempDatCells(dsdx);
  dg_dat_pool->releaseTempDatCells(dsdy);

  for(auto const &p : pointMap) {
    LSKinkPoint point = p.second;
    point.s = point.s / (DG_FP)point.count;
    point.dsdx = point.dsdx / (DG_FP)point.count;
    point.dsdy = point.dsdy / (DG_FP)point.count;
    point.count = 1;
  }

  DG_FP *kink_ptr = getOP2PtrHost(kink, OP_RW);
  for(auto const &p : pointMap) {
    auto const neighbouring_points = point_map_for_kink_detection.find(p.first);
    if(neighbouring_points == point_map_for_kink_detection.end())
      throw std::runtime_error("Error when looking up neighbouring points for kink detection");

    const std::vector<DGUtils::Vec<2>> &points = neighbouring_points->second;

    LSKinkPoint main_point = p.second;
    const DG_FP n1_mag = sqrt(main_point.dsdx * main_point.dsdx + main_point.dsdy * main_point.dsdy);
    const DG_FP n1_x = main_point.dsdx / n1_mag;
    const DG_FP n1_y = main_point.dsdy / n1_mag;
    
    bool is_a_kink = false;
    for(int j = 0; j < points.size(); j++) {
      if(pointMap.count(points[j]) == 0)
        continue;
      LSKinkPoint &point_j = pointMap[points[j]];
      const DG_FP n2_mag = sqrt(point_j.dsdx * point_j.dsdx + point_j.dsdy * point_j.dsdy);
      const DG_FP n2_x = point_j.dsdx / n2_mag;
      const DG_FP n2_y = point_j.dsdy / n2_mag;

      const DG_FP sqr_dist = (n1_x - n2_x) * (n1_x - n2_x) + (n1_y - n2_y) * (n1_y - n2_y);
      if(sqr_dist > kink_sqr_tol) {
        is_a_kink = true;
        break;
      }
    }
    if(is_a_kink) {
      for(auto const &ind : p.second.indices) {
        kink_ptr[ind] = 1.0;
      }
    }
  }
  releaseOP2PtrHost(kink, OP_RW, kink_ptr);
}

#ifdef PRINT_SAMPLE_PTS
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

std::string ls_double_to_text(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}
#endif

void LevelSetSolver2D::reinitLS() {
  timer->startTimer("LevelSetSolver2D - reinitLS");

  timer->startTimer("LevelSetSolver2D - get cells containing interface");
  const DG_FP *s_ptr = getOP2PtrHostHE(s, OP_READ);
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
  timer->endTimer("LevelSetSolver2D - get cells containing interface");

  timer->startTimer("LevelSetSolver2D - create polynomials");
  const DG_FP *_x_ptr = getOP2PtrHostHE(mesh->x, OP_READ);
  const DG_FP *_y_ptr = getOP2PtrHostHE(mesh->y, OP_READ);

  std::map<int,int> _cell2polyMap;
  std::vector<PolyApprox> _polys;
  std::map<int,std::set<int>> stencils = PolyApprox::get_stencils(cellInds, mesh->face2cells, _x_ptr, _y_ptr);

  // Populate map
  int i = 0;
  std::set<int> stencilInds;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    std::set<int> stencil = stencils.at(*it);
    stencilInds.insert(*it);
    stencilInds.insert(stencil.begin(), stencil.end());
    PolyApprox p(*it, stencil, _x_ptr, _y_ptr, s_ptr, h);
    _polys.push_back(p);
    _cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostHE(mesh->x, OP_READ, _x_ptr);
  releaseOP2PtrHostHE(mesh->y, OP_READ, _y_ptr);
  releaseOP2PtrHostHE(s, OP_READ, s_ptr);
  timer->endTimer("LevelSetSolver2D - create polynomials");

  if(kink_detection) {
    timer->startTimer("LevelSetSolver2D - detect kinks");
    detect_kinks(stencilInds);
    timer->endTimer("LevelSetSolver2D - detect kinks");
  }

  DGTempDat tmp_sampleX = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleY = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  timer->startTimer("LevelSetSolver2D - sample interface");
  sampleInterface(tmp_sampleX.dat, tmp_sampleY.dat, _polys, _cell2polyMap, cellInds);
  timer->endTimer("LevelSetSolver2D - sample interface");

  timer->startTimer("LevelSetSolver2D - build KD-Tree");
  const DG_FP *sample_pts_x = getOP2PtrHost(tmp_sampleX.dat, OP_READ);
  const DG_FP *sample_pts_y = getOP2PtrHost(tmp_sampleY.dat, OP_READ);

#ifdef PRINT_SAMPLE_PTS
  std::ofstream file("points.txt");
  file << "X,Y" << std::endl;

  for(int i = 0; i < LS_SAMPLE_NP * mesh->cells->size; i++) {
    if(!isnan(sample_pts_x[i]) && !isnan(sample_pts_y[i])) {
      file << ls_double_to_text(sample_pts_x[i]) << ",";
      file << ls_double_to_text(sample_pts_y[i]) << std::endl;
    }
  }

  file.close();
#endif

  kdtree->reset();
  kdtree->set_poly_data(_polys, _cell2polyMap);
  kdtree->build_tree(sample_pts_x, sample_pts_y, LS_SAMPLE_NP * mesh->cells->size, s);

  releaseOP2PtrHost(tmp_sampleX.dat, OP_READ, sample_pts_x);
  releaseOP2PtrHost(tmp_sampleY.dat, OP_READ, sample_pts_y);

  dg_dat_pool->releaseTempDatCells(tmp_sampleX);
  dg_dat_pool->releaseTempDatCells(tmp_sampleY);
  timer->endTimer("LevelSetSolver2D - build KD-Tree");

  timer->startTimer("LevelSetSolver2D - query KD-Tree");
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);

  DG_FP *closest_x = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_y = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  int *poly_ind    = (int *)calloc(DG_NP * mesh->cells->size, sizeof(int));
  std::vector<PolyApprox> polys;

  DG_FP *surface_ptr = getOP2PtrHost(s, OP_RW);
  if(!kdtree->empty) {
    #pragma omp parallel for
    for(int i = 0; i < mesh->cells->size; i++) {
      // Get closest sample point
      KDCoord tmp = kdtree->closest_point(x_ptr[i * DG_NP], y_ptr[i * DG_NP]);
      const DG_FP dist2 = (tmp.coord[0] - x_ptr[i * DG_NP]) * (tmp.coord[0] - x_ptr[i * DG_NP])
                        + (tmp.coord[1] - y_ptr[i * DG_NP]) * (tmp.coord[1] - y_ptr[i * DG_NP]);
      closest_x[i * DG_NP] = tmp.coord[0];
      closest_y[i * DG_NP] = tmp.coord[1];
      poly_ind[i * DG_NP]  = tmp.poly;
      if(dist2 < (reinit_width + 2.0 * h) * (reinit_width + 2.0 * h)) {
        for(int j = 1; j < DG_NP; j++) {
          // Get closest sample point
          KDCoord tmp = kdtree->closest_point(x_ptr[i * DG_NP + j], y_ptr[i * DG_NP + j]);
          closest_x[i * DG_NP + j] = tmp.coord[0];
          closest_y[i * DG_NP + j] = tmp.coord[1];
          poly_ind[i * DG_NP + j]  = tmp.poly;
        }
      }
    }

    polys = kdtree->get_polys();
  } else {
    #pragma omp parallel for
    for(int i = 0; i < mesh->cells->size * DG_NP; i++) {
      if(surface_ptr[i] > 0.0) {
        surface_ptr[i] = ls_cap;
      } else {
        surface_ptr[i] = -ls_cap;
      }
    }
  }
  timer->endTimer("LevelSetSolver2D - query KD-Tree");

  timer->startTimer("LevelSetSolver2D - newton method");
  // Newton method
  DG_FP *kink_ptr = getOP2PtrHost(kink, OP_RW);
  if(!kdtree->empty) {
    newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, x_ptr,
                  y_ptr, poly_ind, polys, surface_ptr, kink_ptr, h, reinit_width, 
                  ls_cap, kink_avoid_whole_element);
  }
  releaseOP2PtrHost(kink, OP_RW, kink_ptr);
  releaseOP2PtrHost(s, OP_RW, surface_ptr);
  timer->endTimer("LevelSetSolver2D - newton method");

  free(closest_x);
  free(closest_y);
  free(poly_ind);

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  timer->endTimer("LevelSetSolver2D - reinitLS");
}
