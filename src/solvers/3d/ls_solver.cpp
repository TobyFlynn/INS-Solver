#include "solvers/3d/ls_solver.h"

#include "op_seq.h"

#include "ls_utils/3d/kd_tree.h"
#ifdef INS_MPI
#include "ls_utils/3d/kd_tree_mpi.h"
#endif
#include "op2_utils.h"
#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "dg_dat_pool.h"

#include <iostream>
#include <fstream>
#include <limits>

#include "timing.h"
#include "config.h"

extern Timing *timer;
extern DGDatPool *dg_dat_pool;
extern Config *config;

/**************************************************************************
 * LS Advection Solver class that extends the base Advection Solver class *
 **************************************************************************/
LevelSetAdvectionSolver3D::LevelSetAdvectionSolver3D(DGMesh3D *m) : AdvectionSolver3D(m) {}

void LevelSetAdvectionSolver3D::set_bc_types(op_dat bc) {
  bc_types = bc;
}

void LevelSetAdvectionSolver3D::bc_kernel(op_dat val, op_dat u, op_dat v,
                                          op_dat w, op_dat out) {
  op_par_loop(ls_advec_3d_bc, "ls_advec_3d_bc", mesh->bfaces,
              op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
}

void LevelSetAdvectionSolver3D::bc_kernel_oi(op_dat val, op_dat u, op_dat v,
                                             op_dat w, op_dat flux) {
  op_par_loop(ls_advec_3d_oi_bc, "ls_advec_3d_oi_bc", mesh->bfaces,
              op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_INC));
}

/************************
 * Main LS Solver class *
 ************************/
LevelSetSolver3D::LevelSetSolver3D(DGMesh3D *m) {
  mesh = m;
  resuming = false;
  advectionSolver = new LevelSetAdvectionSolver3D(m);

  s = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_solver_s");
  kink = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_solver_kink");

  h = 0;
}

LevelSetSolver3D::LevelSetSolver3D(DGMesh3D *m, const std::string &filename) {
  mesh = m;
  resuming = true;
  advectionSolver = new LevelSetAdvectionSolver3D(m);

  s = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ls_solver_s");
  kink = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_solver_kink");

  h = 0;
}

LevelSetSolver3D::~LevelSetSolver3D() {
  delete advectionSolver;
  delete kdtree;
}

void LevelSetSolver3D::set_bc_types(op_dat bc) {
  advectionSolver->set_bc_types(bc);
}

void LevelSetSolver3D::init() {
  if(!resuming) {
    op_par_loop(init_surface_3d, "init_surface_3d", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  h = std::numeric_limits<DG_FP>::max();
  op_par_loop(calc_h_3d_max, "calc_h_3d_max", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MIN));
  h = 1.0 / h;

/*
  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
*/

  op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
              op_arg_dat(kink, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // alpha = 2.0 * h / DG_ORDER;
  // order_width = 2.0 * h;
  // epsilon = h / DG_ORDER;
  alpha = 24.0 * h;
  // order_width = 12.0 * h;
  // epsilon = h;
  ls_cap = 100.0;
  // reinit_width = fmin(4.0 * alpha, ls_cap);
  reinit_width = ls_cap;
  // reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  // numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  reinit_counter = 0;
  reinit_frequency = 24;
  config->getInt("level-set-options", "reinit_freq", reinit_frequency);
  reinitialise = reinit_frequency > 0;

  op_printf("LS h: %g\nLS alpha: %g\n", h, alpha);

  #ifdef INS_MPI
  kdtree = new KDTree3DMPI(mesh, 2.0 * reinit_width);
  #else
  kdtree = new KDTree3D(mesh);
  #endif

  int tmp_reinit_ic = 1;
  config->getInt("level-set-options", "reinit_ic", tmp_reinit_ic);

  int tmp_kink = 1;
  config->getInt("level-set-options", "kink_detection", tmp_kink);
  kink_detection = tmp_kink == 1;

  if(tmp_reinit_ic != 0)
    reinitLS();

  op_par_loop(ls_post_reinit, "ls_post_reinit", mesh->cells,
                op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

}

void LevelSetSolver3D::getRhoMu(op_dat rho, op_dat mu) {
  timer->startTimer("LevelSetSolver3D - getRhoMu");
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver3D - getRhoMu");
}

void LevelSetSolver3D::getRhoVolOI(op_dat rho) {
  timer->startTimer("LevelSetSolver3D - getRhoVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, s, 0.0, rho);
  op_par_loop(ls_3d_step_rho_vol_oi, "ls_3d_step_rho_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver3D - getRhoVolOI");
}

void LevelSetSolver3D::getMuVolOI(op_dat mu) {
  timer->startTimer("LevelSetSolver3D - getMuVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, s, 0.0, mu);
  op_par_loop(ls_3d_step_mu_vol_oi, "ls_3d_step_mu_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(mu, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver3D - getMuVolOI");
}

void LevelSetSolver3D::getRhoSurfOI(op_dat rho) {
  timer->startTimer("LevelSetSolver3D - getRhoSurfOI");
  DGTempDat rho_tmp = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  op_par_loop(ls_step_surf_oi_0, "ls_step_surf_oi_0", mesh->cells,
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho_tmp.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF3D_INTERP, rho_tmp.dat, 0.0, rho);
  dg_dat_pool->releaseTempDatCells(rho_tmp);
  op_par_loop(ls_3d_step_surf_oi_1, "ls_3d_step_surf_oi_1", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver3D - getRhoSurfOI");
}

void LevelSetSolver3D::getNormalsCurvature(op_dat nx, op_dat ny, op_dat nz, op_dat curv) {
  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  timer->startTimer("LevelSetSolver3D - getNormalsCurvature");
  mesh->grad(s, nx, ny, nz);
  mesh->div(nx, ny, nz, curv);
  timer->endTimer("LevelSetSolver3D - getNormalsCurvature");
}

void LevelSetSolver3D::step(op_dat u, op_dat v, op_dat w, const DG_FP dt, const int num_steps) {
  timer->startTimer("LevelSetSolver3D - step");
  advectionSolver->set_dt(dt);
  for(int i = 0; i < num_steps; i++)
    advectionSolver->step(s, u, v, w);

  reinit_counter++;
  if(reinitialise && reinit_counter >= reinit_frequency) {
    reinitLS();
    reinit_counter = 0;
    op_par_loop(ls_post_reinit, "ls_post_reinit", mesh->cells,
                op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }
  timer->endTimer("LevelSetSolver3D - step");
}

bool newtoncp_gepp(arma::mat &A, arma::vec &b) {
  for (int i = 0; i < 4; ++i) {
    int j = i;
    for (int k = i + 1; k < 4; ++k)
      if (std::abs(A(k,i)) > std::abs(A(j,i)))
        j = k;
    if (j != i) {
      for (int k = 0; k < 4; ++k)
        std::swap(A(i,k), A(j,k));
      std::swap(b(i), b(j));
    }

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<DG_FP>::epsilon())
      return false;

    DG_FP fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 4; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 4; ++j) {
      for (int k = i + 1; k < 4; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 4 - 1; i >= 0; --i) {
    DG_FP sum = 0.0;
    for (int j = i + 1; j < 4; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(DG_FP &closest_pt_x, DG_FP &closest_pt_y, DG_FP &closest_pt_z,
                   const DG_FP node_x, const DG_FP node_y, const DG_FP node_z,
                   PolyApprox3D &p, const DG_FP h) {
  DG_FP lambda = 0.0;
  bool converged = false;
  DG_FP pt_x = closest_pt_x;
  DG_FP pt_y = closest_pt_y;
  DG_FP pt_z = closest_pt_z;
  DG_FP init_x = closest_pt_x;
  DG_FP init_y = closest_pt_y;
  DG_FP init_z = closest_pt_z;

  for(int step = 0; step < 50; step++) {
    DG_FP pt_x_old = pt_x;
    DG_FP pt_y_old = pt_y;
    DG_FP pt_z_old = pt_z;
    // Evaluate surface and gradient at current guess
    DG_FP surface = p.val_at(pt_x, pt_y, pt_z);
    DG_FP surface_dx, surface_dy, surface_dz;
    p.grad_at(pt_x, pt_y, pt_z, surface_dx, surface_dy, surface_dz);
    // Evaluate Hessian
    DG_FP hessian[6];
    p.hessian_at(pt_x, pt_y, pt_z, hessian[0], hessian[1], hessian[2],
                 hessian[3], hessian[4], hessian[5]);

    // Check if |nabla(surface)| = 0, if so then return
    DG_FP gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy + surface_dz * surface_dz;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy + (node_z - pt_z) * surface_dz) / gradsqrnorm;

    // Gradient of functional
    arma::vec gradf(4);
    gradf(0) = pt_x - node_x + lambda * surface_dx;
    gradf(1) = pt_y - node_y + lambda * surface_dy;
    gradf(2) = pt_z - node_z + lambda * surface_dz;
    gradf(3) = surface;

    // Calculate Hessian of functional
    arma::mat hessianf(4, 4);
    hessianf(0, 0) = 1.0 + lambda * hessian[0];
    hessianf(0, 1) = lambda * hessian[3]; hessianf(1, 0) = hessianf(0, 1);
    hessianf(0, 2) = lambda * hessian[4]; hessianf(2, 0) = hessianf(0, 2);
    hessianf(0, 3) = surface_dx; hessianf(3, 0) = hessianf(0, 3);

    hessianf(1, 1) = 1.0 + lambda * hessian[1];
    hessianf(1, 2) = lambda * hessian[5]; hessianf(2, 1) = hessianf(1, 2);
    hessianf(1, 3) = surface_dy; hessianf(3, 1) = hessianf(1, 3);

    hessianf(2, 2) = 1.0 + lambda * hessian[2];
    hessianf(2, 3) = surface_dz; hessianf(3, 2) = hessianf(2, 3);

    hessianf(3, 3) = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      DG_FP delta1_x = (surface / gradsqrnorm) * surface_dx;
      DG_FP delta1_y = (surface / gradsqrnorm) * surface_dy;
      DG_FP delta1_z = (surface / gradsqrnorm) * surface_dz;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy + (node_z - pt_z) * surface_dz) / gradsqrnorm;
      DG_FP delta2_x = pt_x - node_x + lambda * surface_dx;
      DG_FP delta2_y = pt_y - node_y + lambda * surface_dy;
      DG_FP delta2_z = pt_z - node_z + lambda * surface_dz;
      DG_FP msqr = delta2_x * delta2_x + delta2_y * delta2_y + delta2_z * delta2_z;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
        delta2_z *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
      pt_z -= delta1_z + delta2_z;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      DG_FP msqr = ans(0) * ans(0) + ans(1) * ans(1) + ans(2) * ans(2);
      if(msqr > h * 0.5 * h * 0.5)
        ans = ans * 0.5 * h / sqrt(msqr);

      // Update guess
      pt_x -= ans(0);
      pt_y -= ans(1);
      pt_z -= ans(2);
      lambda -= ans(3);
    }

    // Gone outside the element, return
    if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) + (init_z - pt_z) * (init_z - pt_z) > 1.5 * 1.5 * h * h) {
      pt_x = pt_x_old;
      pt_y = pt_y_old;
      pt_z = pt_z_old;
      break;
    }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) + (pt_z_old - pt_z) * (pt_z_old - pt_z) < 1e-24) {
      converged = true;
      break;
    }
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;
  closest_pt_z = pt_z;

  return converged;
}

void newton_method(const int numPts, DG_FP *closest_x, DG_FP *closest_y, DG_FP *closest_z,
                   const DG_FP *x, const DG_FP *y, const DG_FP *z, int *poly_ind,
                   std::vector<PolyApprox3D> &polys, DG_FP *s, DG_FP *kink, const DG_FP h,
                   const DG_FP reinit_width, const DG_FP ls_cap) {
  int numNonConv = 0;
  int numReinit = 0;

  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    DG_FP dist2 = (closest_x[start_ind] - x[start_ind]) * (closest_x[start_ind] - x[start_ind])
                + (closest_y[start_ind] - y[start_ind]) * (closest_y[start_ind] - y[start_ind])
                + (closest_z[start_ind] - z[start_ind]) * (closest_z[start_ind] - z[start_ind]);
    if(dist2 < (reinit_width + 2.0 * h) * (reinit_width + 2.0 * h)) {
      DG_FP dist1 = (closest_x[i] - x[i]) * (closest_x[i] - x[i])
                  + (closest_y[i] - y[i]) * (closest_y[i] - y[i])
                  + (closest_z[i] - z[i]) * (closest_z[i] - z[i]);
      if(kink[i] != 0.0 && dist1 < 9.0 * h * h)
        continue;
      else
        kink[i] = 0.0;
      
      DG_FP off_x, off_y, off_z;
      polys[poly_ind[i]].get_offsets(off_x, off_y, off_z);
      DG_FP _closest_x = closest_x[i] - off_x;
      DG_FP _closest_y = closest_y[i] - off_y;
      DG_FP _closest_z = closest_z[i] - off_z;
      DG_FP _node_x = x[i] - off_x;
      DG_FP _node_y = y[i] - off_y;
      DG_FP _node_z = z[i] - off_z;
      bool converged = newton_kernel(_closest_x, _closest_y, _closest_z, _node_x, _node_y, _node_z, polys[poly_ind[i]], h);
      if(converged) {
        bool negative = s[i] < 0.0;
        s[i] = (_closest_x - _node_x) * (_closest_x - _node_x) + (_closest_y - _node_y) * (_closest_y - _node_y) + (_closest_z - _node_z) * (_closest_z - _node_z);
        s[i] = sqrt(s[i]);
        if(negative) s[i] *= -1.0;
      } else {
        bool negative = s[i] < 0.0;
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
  DGUtils::Vec<3> coord;
  DG_FP s, dsdx, dsdy, dsdz;
  int count;
  std::vector<int> indices;
};

void LevelSetSolver3D::detect_kinks() {
  op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
              op_arg_dat(kink, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat dsdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dsdy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dsdz = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(s, dsdx.dat, dsdy.dat, dsdz.dat);

  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHost(mesh->z, OP_READ);
  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  const DG_FP *dsdx_ptr = getOP2PtrHost(dsdx.dat, OP_READ);
  const DG_FP *dsdy_ptr = getOP2PtrHost(dsdy.dat, OP_READ);
  const DG_FP *dsdz_ptr = getOP2PtrHost(dsdz.dat, OP_READ);

  std::map<DGUtils::Vec<3>, LSKinkPoint> pointMap;

  for(int i = 0; i < mesh->cells->size * DG_NP; i++) {
    if(fabs(s_ptr[i]) < 4.0 * h) {
      DGUtils::Vec<3> coord(x_ptr[i], y_ptr[i], z_ptr[i]);
      LSKinkPoint point;
      auto res = pointMap.insert({coord, point});
      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.s = s_ptr[i];
        res.first->second.dsdx = dsdx_ptr[i];
        res.first->second.dsdy = dsdy_ptr[i];
        res.first->second.dsdz = dsdz_ptr[i];
        res.first->second.count = 1;
        res.first->second.indices.push_back(i);
      } else {
        // Point already exists
        res.first->second.s += s_ptr[i];
        res.first->second.dsdx += dsdx_ptr[i];
        res.first->second.dsdy += dsdy_ptr[i];
        res.first->second.dsdz += dsdz_ptr[i];
        res.first->second.count++;
        res.first->second.indices.push_back(i);
      }
    }
  }

  releaseOP2PtrHost(dsdz.dat, OP_READ, dsdz_ptr);
  releaseOP2PtrHost(dsdy.dat, OP_READ, dsdy_ptr);
  releaseOP2PtrHost(dsdx.dat, OP_READ, dsdx_ptr);
  releaseOP2PtrHost(s, OP_READ, s_ptr);
  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);

  dg_dat_pool->releaseTempDatCells(dsdx);
  dg_dat_pool->releaseTempDatCells(dsdy);
  dg_dat_pool->releaseTempDatCells(dsdz);

  // Create a map consisting of buckets that contain a vector of points.
  // Each bucket is a square of length max_distance_between_points.
  // Key to the bucket is the coordinates of its bottom left vertex.
  // const DG_FP max_distance_between_points = 2.0 * h;
  const DG_FP max_distance_between_points = h;
  const DG_FP max_distance_between_points_sqr = max_distance_between_points * max_distance_between_points;
  std::map<DGUtils::Vec<3>, std::vector<LSKinkPoint>> bucket_map;
  for(auto const &p : pointMap) {
    LSKinkPoint point = p.second;
    point.s = point.s / (DG_FP)point.count;
    point.dsdx = point.dsdx / (DG_FP)point.count;
    point.dsdy = point.dsdy / (DG_FP)point.count;
    point.dsdz = point.dsdz / (DG_FP)point.count;
    point.count = 1;

    DGUtils::Vec<3> bucket_coord = point.coord;
    bucket_coord[0] = static_cast<int>(bucket_coord[0] / max_distance_between_points) * max_distance_between_points;
    bucket_coord[1] = static_cast<int>(bucket_coord[1] / max_distance_between_points) * max_distance_between_points;
    bucket_coord[2] = static_cast<int>(bucket_coord[2] / max_distance_between_points) * max_distance_between_points;

    auto res = bucket_map.insert({bucket_coord, std::vector<LSKinkPoint>()});
    res.first->second.push_back(point);
  }

  DGUtils::Vec<3> bucket_offsets[26] = {
    DGUtils::Vec<3>(-max_distance_between_points, 0.0, 0.0),
    DGUtils::Vec<3>(-max_distance_between_points, max_distance_between_points, 0.0),
    DGUtils::Vec<3>(0.0, max_distance_between_points, 0.0),
    DGUtils::Vec<3>(max_distance_between_points, max_distance_between_points, 0.0),
    DGUtils::Vec<3>(max_distance_between_points, 0.0, 0.0),
    DGUtils::Vec<3>(max_distance_between_points, -max_distance_between_points, 0.0),
    DGUtils::Vec<3>(0.0, -max_distance_between_points, 0.0),
    DGUtils::Vec<3>(-max_distance_between_points, -max_distance_between_points, 0.0),
    //
    DGUtils::Vec<3>(-max_distance_between_points, 0.0, max_distance_between_points),
    DGUtils::Vec<3>(-max_distance_between_points, max_distance_between_points, max_distance_between_points),
    DGUtils::Vec<3>(0.0, max_distance_between_points, 0.0),
    DGUtils::Vec<3>(max_distance_between_points, max_distance_between_points, max_distance_between_points),
    DGUtils::Vec<3>(max_distance_between_points, 0.0, 0.0),
    DGUtils::Vec<3>(max_distance_between_points, -max_distance_between_points, max_distance_between_points),
    DGUtils::Vec<3>(0.0, -max_distance_between_points, 0.0),
    DGUtils::Vec<3>(-max_distance_between_points, -max_distance_between_points, max_distance_between_points),
    DGUtils::Vec<3>(0.0, 0.0, max_distance_between_points),
    //
    DGUtils::Vec<3>(-max_distance_between_points, 0.0, -max_distance_between_points),
    DGUtils::Vec<3>(-max_distance_between_points, max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(0.0, max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(max_distance_between_points, max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(max_distance_between_points, 0.0, -max_distance_between_points),
    DGUtils::Vec<3>(max_distance_between_points, -max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(0.0, -max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(-max_distance_between_points, -max_distance_between_points, -max_distance_between_points),
    DGUtils::Vec<3>(0.0, 0.0, -max_distance_between_points)
  };

  DG_FP *kink_ptr = getOP2PtrHost(kink, OP_RW);
  const DG_FP sqr_tol = 0.5 * 0.5;
  for(auto const &bucket : bucket_map) {
    // Get points in this bucket
    std::vector<LSKinkPoint> points(bucket.second);
    // Get points in surrounding buckets
    for(int i = 0; i < 26; i++) {
      DGUtils::Vec<3> bucket_coord = bucket.first + bucket_offsets[i];
      auto const neighbouring_bucket = bucket_map.find(bucket_coord);
      if(neighbouring_bucket != bucket_map.end()) {
        points.insert(points.end(), neighbouring_bucket->second.begin(), neighbouring_bucket->second.end());
      }
    }

    // Iterate over all points
    for(int i = 0; i < points.size(); i++) {
      bool is_a_kink = false;
      for(int j = 0; j < points.size(); j++) {
        if(j == i || (points[i].coord - points[j].coord).sqr_magnitude() > max_distance_between_points_sqr)
          continue;

        const DG_FP n1_mag = sqrt(points[j].dsdx * points[j].dsdx + points[j].dsdy * points[j].dsdy + points[j].dsdz * points[j].dsdz);
        const DG_FP n1_x = points[j].dsdx / n1_mag;
        const DG_FP n1_y = points[j].dsdy / n1_mag;
        const DG_FP n1_z = points[j].dsdz / n1_mag;
        for(int k = j + 1; k < points.size(); k++) {
          if(k == i || k == j || (points[i].coord - points[k].coord).sqr_magnitude() > max_distance_between_points_sqr)
            continue;

          const DG_FP n2_mag = sqrt(points[k].dsdx * points[k].dsdx + points[k].dsdy * points[k].dsdy + points[k].dsdz * points[k].dsdz);
          const DG_FP n2_x = points[k].dsdx / n2_mag;
          const DG_FP n2_y = points[k].dsdy / n2_mag;
          const DG_FP n2_z = points[k].dsdz / n2_mag;

          const DG_FP sqr_dist = (n1_x - n2_x) * (n1_x - n2_x) + (n1_y - n2_y) * (n1_y - n2_y) + (n1_z - n2_z) * (n1_z - n2_z);
          if(sqr_dist > sqr_tol) {
            is_a_kink = true;
            break;
          }
        }
        if(is_a_kink)
          break;
      }

      if(is_a_kink) {
        for(auto const &ind : points[i].indices) {
          kink_ptr[ind] = 1.0;
        }
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

void LevelSetSolver3D::reinitLS() {
  timer->startTimer("LevelSetSolver3D - reinitLS");
  if(kink_detection) {
    timer->startTimer("LevelSetSolver2D - detect kinks");
    detect_kinks();
    timer->endTimer("LevelSetSolver2D - detect kinks");
  }

  if(h == 0.0) {
    op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
                op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
                op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
    h = 1.0 / h;
  }

  timer->startTimer("LevelSetSolver3D - get cells containing interface");
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
  timer->endTimer("LevelSetSolver3D - get cells containing interface");

  timer->startTimer("LevelSetSolver3D - create polynomials");
  std::map<int,int> _cell2polyMap;
  std::vector<PolyApprox3D> _polys;

  const DG_FP *_x_ptr = getOP2PtrHostHE(mesh->x, OP_READ);
  const DG_FP *_y_ptr = getOP2PtrHostHE(mesh->y, OP_READ);
  const DG_FP *_z_ptr = getOP2PtrHostHE(mesh->z, OP_READ);

  std::map<int,std::set<int>> stencils = PolyApprox3D::get_stencils(cellInds, mesh->face2cells, _x_ptr, _y_ptr, _z_ptr);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    std::set<int> stencil = stencils.at(*it);
    PolyApprox3D p(*it, stencil, _x_ptr, _y_ptr, _z_ptr, s_ptr, h);
    _polys.push_back(p);
    _cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostHE(mesh->x, OP_READ, _x_ptr);
  releaseOP2PtrHostHE(mesh->y, OP_READ, _y_ptr);
  releaseOP2PtrHostHE(mesh->z, OP_READ, _z_ptr);
  releaseOP2PtrHostHE(s, OP_READ, s_ptr);
  timer->endTimer("LevelSetSolver3D - create polynomials");

  DGTempDat tmp_sampleX = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleY = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleZ = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  timer->startTimer("LevelSetSolver3D - sample interface");
  sampleInterface(tmp_sampleX.dat, tmp_sampleY.dat, tmp_sampleZ.dat, _polys, _cell2polyMap, cellInds);
  timer->endTimer("LevelSetSolver3D - sample interface");

  timer->startTimer("LevelSetSolver3D - build KD-Tree");
  const DG_FP *sample_pts_x = getOP2PtrHost(tmp_sampleX.dat, OP_READ);
  const DG_FP *sample_pts_y = getOP2PtrHost(tmp_sampleY.dat, OP_READ);
  const DG_FP *sample_pts_z = getOP2PtrHost(tmp_sampleZ.dat, OP_READ);

#ifdef PRINT_SAMPLE_PTS
  std::ofstream file("points.txt");
  file << "X,Y,Z" << std::endl;

  for(int i = 0; i < LS_SAMPLE_NP * mesh->cells->size; i++) {
    if(!isnan(sample_pts_x[i]) && !isnan(sample_pts_y[i]) && !isnan(sample_pts_z[i])) {
      file << ls_double_to_text(sample_pts_x[i]) << ",";
      file << ls_double_to_text(sample_pts_y[i]) << ",";
      file << ls_double_to_text(sample_pts_z[i]) << std::endl;
    }
  }

  file.close();
#endif

  kdtree->reset();
  kdtree->set_poly_data(_polys, _cell2polyMap);
  kdtree->build_tree(sample_pts_x, sample_pts_y, sample_pts_z, LS_SAMPLE_NP * mesh->cells->size, s);

  releaseOP2PtrHost(tmp_sampleX.dat, OP_READ, sample_pts_x);
  releaseOP2PtrHost(tmp_sampleY.dat, OP_READ, sample_pts_y);
  releaseOP2PtrHost(tmp_sampleZ.dat, OP_READ, sample_pts_z);

  dg_dat_pool->releaseTempDatCells(tmp_sampleX);
  dg_dat_pool->releaseTempDatCells(tmp_sampleY);
  dg_dat_pool->releaseTempDatCells(tmp_sampleZ);
  timer->endTimer("LevelSetSolver3D - build KD-Tree");

  timer->startTimer("LevelSetSolver3D - query KD-Tree");
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHost(mesh->z, OP_READ);

  DG_FP *closest_x = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_y = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_z = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  int *poly_ind     = (int *)calloc(DG_NP * mesh->cells->size, sizeof(int));
  std::vector<PolyApprox3D> polys;

  DG_FP *surface_ptr = getOP2PtrHost(s, OP_RW);
  if(!kdtree->empty) {
    #pragma omp parallel for
    for(int i = 0; i < mesh->cells->size; i++) {
      // Get closest sample point
      KDCoord tmp = kdtree->closest_point(x_ptr[i * DG_NP], y_ptr[i * DG_NP], z_ptr[i * DG_NP]);
      const DG_FP dist2 = (tmp.coord[0] - x_ptr[i * DG_NP]) * (tmp.coord[0] - x_ptr[i * DG_NP])
                        + (tmp.coord[1] - y_ptr[i * DG_NP]) * (tmp.coord[1] - y_ptr[i * DG_NP])
                        + (tmp.coord[2] - z_ptr[i * DG_NP]) * (tmp.coord[2] - z_ptr[i * DG_NP]);
      closest_x[i * DG_NP] = tmp.coord[0];
      closest_y[i * DG_NP] = tmp.coord[1];
      closest_z[i * DG_NP] = tmp.coord[2];
      poly_ind[i * DG_NP]  = tmp.poly;
      if(dist2 < (reinit_width + 2.0 * h) * (reinit_width + 2.0 * h)) {
        for(int j = 1; j < DG_NP; j++) {
          // Get closest sample point
          KDCoord tmp = kdtree->closest_point(x_ptr[i * DG_NP + j], y_ptr[i * DG_NP + j], z_ptr[i * DG_NP + j]);
          closest_x[i * DG_NP + j] = tmp.coord[0];
          closest_y[i * DG_NP + j] = tmp.coord[1];
          closest_z[i * DG_NP + j] = tmp.coord[2];
          poly_ind[i * DG_NP + j]  = tmp.poly;
        }
      }
    }

    polys = kdtree->get_polys();
  }
  timer->endTimer("LevelSetSolver3D - query KD-Tree");

  timer->startTimer("LevelSetSolver3D - newton method");
  // Newton method
  DG_FP *kink_ptr = getOP2PtrHost(kink, OP_RW);
  if(!kdtree->empty) {
    newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, closest_z,
                  x_ptr, y_ptr, z_ptr, poly_ind, polys, surface_ptr, kink_ptr, 
                  h, reinit_width, ls_cap);
  }
  releaseOP2PtrHost(kink, OP_RW, kink_ptr);
  releaseOP2PtrHost(s, OP_RW, surface_ptr);
  timer->endTimer("LevelSetSolver3D - newton method");

  free(closest_x);
  free(closest_y);
  free(closest_z);
  free(poly_ind);

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);
  timer->endTimer("LevelSetSolver3D - reinitLS");
}
