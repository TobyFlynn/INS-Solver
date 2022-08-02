#include "ls.h"

#include "op_seq.h"

#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>

#include "dg_constants.h"
#include "dg_utils.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "dg_compiler_defs.h"

#include "kd_tree.h"
#include "timing.h"
#include "ls_reinit_stencil.h"

extern Timing *timer;

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

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<double>::epsilon())
      return false;

    double fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 3; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 3; ++j) {
      for (int k = i + 1; k < 3; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 3 - 1; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < 3; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(double &closest_pt_x, double &closest_pt_y,
                   const double node_x, const double node_y,
                   const int cell_ind, op_map edge_map, op_dat x_dat, op_dat y_dat,
                   op_dat s_dat, const double h) {
  double lambda = 0.0;
  bool converged = false;
  double pt_x = closest_pt_x;
  double pt_y = closest_pt_y;
  double init_x = closest_pt_x;
  double init_y = closest_pt_y;

  PolyApprox p(3, cell_ind, edge_map, x_dat, y_dat, s_dat);

  for(int step = 0; step < 10; step++) {
    double pt_x_old = pt_x;
    double pt_y_old = pt_y;
    // Evaluate surface and gradient at current guess
    double surface = p.val_at(pt_x, pt_y);
    double surface_dx, surface_dy;
    p.grad_at(pt_x, pt_y, surface_dx, surface_dy);
    // Evaluate Hessian
    double hessian[3];
    p.hessian_at(pt_x, pt_y, hessian[0], hessian[1], hessian[2]);
    
    // Check if |nabla(surface)| = 0, if so then return
    double gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy;
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
      double delta1_x = (surface / gradsqrnorm) * surface_dx;
      double delta1_y = (surface / gradsqrnorm) * surface_dy;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy) / gradsqrnorm;
      double delta2_x = pt_x - node_x + lambda * surface_dx;
      double delta2_y = pt_y - node_y + lambda * surface_dy;
      double msqr = delta2_x * delta2_x + delta2_y * delta2_y;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      double msqr = ans(0) * ans(0) + ans(1) * ans(1);
      if(msqr > h * 0.5 * h * 0.5)
        ans = ans * 0.5 * h / sqrt(msqr);

      // Update guess
      pt_x -= ans(0);
      pt_y -= ans(1);
      lambda -= ans(2);
    }

    // Gone outside the element, return
    // if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
    //   pt_x = pt_x_old;
    //   pt_y = pt_y_old;
    //   break;
    // }

    // Converged, no more steps required
    if((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) < 1e-18) {
      converged = true;
      break;
    }
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;

  return converged;
}

void newton_method(const int numPts, double *closest_x, double *closest_y, const double *x, const double *y,
                   int *cell_ind, op_map edge_map, op_dat x_dat, op_dat y_dat, op_dat s_dat, const double h) {
  int numNonConv = 0;
  int numReinit = 0;
  double *s_ptr = (double *)s_dat->data;
  double *s_new = (double *)calloc(numPts, sizeof(double));
  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s_ptr[start_ind + j]) < 0.05) {
        reinit = true;
      }
    }

    s_new[i] = s_ptr[i];

    if(reinit) {
      bool tmp = newton_kernel(closest_x[i], closest_y[i], x[i], y[i], cell_ind[i], 
                               edge_map, x_dat, y_dat, s_dat, h);
      if(tmp) {
        bool negative = s_ptr[i] < 0.0;
        s_new[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]);
        s_new[i] = sqrt(s_new[i]);
        if(negative) s_new[i] *= -1.0;
      }
      if(!tmp) {
        #pragma omp atomic
        numNonConv++;
      }
      #pragma omp atomic
      numReinit++;
    }
  }

  memcpy(s_ptr, s_new, numPts * sizeof(double));

  free(s_new);
  
  if(numNonConv != 0 || numReinit == 0)
    std::cout << numNonConv << " non-converged points out of " << numReinit << " points reinitialised" << std::endl;
}

void LS::reinit_ls() {
  timer->startTimer("LS - Sample Interface");
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s, 0.0, s_modal);

  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s_sample_x,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE),
              op_arg_dat(s_sample_y,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE));

  timer->endTimer("LS - Sample Interface");

  op_arg op2_args[] = {
    op_arg_dat(s_sample_x, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(s_sample_y, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(s, -1, OP_ID, DG_NP, "double", OP_RW),
    op_arg_dat(s_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges(s_sample_x->set, 12, op2_args);

  std::string fileName = "interface_sample_points.txt";
  std::ofstream out_file(fileName.c_str());
  out_file << "X,Y,Z" << std::endl;

  const double *x_coords = (double *)s_sample_x->data;
  const double *y_coords = (double *)s_sample_y->data;
  for(int i = 0; i < LS_SAMPLE_NP * mesh->numCells; i++) {
    if(!isnan(x_coords[i]) && !isnan(y_coords[i]))
      out_file << x_coords[i] << "," << y_coords[i] << ",0.0" << std::endl;
  }

  out_file.close();

  timer->startTimer("LS - Construct K-D Tree");
  KDTree kdtree((double *)s_sample_x->data, (double *)s_sample_y->data, LS_SAMPLE_NP * mesh->numCells);
  timer->endTimer("LS - Construct K-D Tree");

  const double *mesh_x_coords = (double *)mesh->x->data;
  const double *mesh_y_coords = (double *)mesh->y->data;

  timer->startTimer("LS - Newton Method");
  double *closest_x = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  double *closest_y = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  int *cell_ind     = (int *)calloc(DG_NP * mesh->numCells, sizeof(int));

  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->numCells; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(mesh_x_coords[i], mesh_y_coords[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    cell_ind[i]  = tmp.cell;
  }

  newton_method(DG_NP * mesh->numCells, closest_x, closest_y, mesh_x_coords, mesh_y_coords,
                cell_ind, mesh->edge2cells, mesh->x, mesh->y, s, h);

  free(closest_x);
  free(closest_y);
  free(cell_ind);
  timer->endTimer("LS - Newton Method");

  op_mpi_set_dirtybit(12, op2_args);
}
