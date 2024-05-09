#include "ls_utils/2d/ls_reinit_poly.h"

#include <random>
#include <map>
#include <fstream>
#include <iostream>

#include "op2_utils.h"
#include "timing.h"
#include "dg_constants/dg_constants_2d.h"
#include "dg_global_constants/dg_global_constants_2d.h"

#include "CDT.h"

extern Timing *timer;
extern DGConstants *constants;

using namespace std;

struct Point {
  DGUtils::Vec<2> coord;
  DG_FP value;
  int count;
};

void get_stencil_values(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                        const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP offset_x,
                        const DG_FP offset_y, map<DGUtils::Vec<2>, Point> &pointMap);

bool vec_contains(const int val, const vector<int> &vec) {
  for(int i = 0; i < vec.size(); i++) {
    if(val == vec[i])
      return true;
  }
  return false;
}

void PolyApprox::calc_offset(const int ind, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  offset_x = 0.0;
  offset_y = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    offset_x += x_ptr[ind * DG_NP + i];
    offset_y += y_ptr[ind * DG_NP + i];
  }
  offset_x /= (DG_FP)DG_NP;
  offset_y /= (DG_FP)DG_NP;
}

void PolyApprox::stencil_data(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                              const DG_FP *y_ptr, const DG_FP *s_ptr, vector<DG_FP> &x,
                              vector<DG_FP> &y, vector<DG_FP> &s_mean) {
  map<DGUtils::Vec<2>, Point> pointMap;
  get_stencil_values(cell_ind, stencil, x_ptr, y_ptr, s_ptr, offset_x, offset_y, pointMap);

  for(auto const &p : pointMap) {
    x.push_back(p.second.coord[0]);
    y.push_back(p.second.coord[1]);
    s_mean.push_back(p.second.value / (DG_FP)p.second.count);
    /*
    // Calculate mean
    DG_FP sum = 0.0;
    for(const auto &val : p.second.values) {
      sum += val;
    }
    const DG_FP mean = sum / (DG_FP)p.second.values.size();
    s_mean.push_back(mean);
    // Calculate variance
    if(p.second.values.size() == 1) {
      s_variance.push_back(0.0);
    } else {
      sum = 0.0;
      for(const auto &val : p.second.values) {
        sum += (val - mean) * (val - mean);
      }
      const DG_FP variance = sum / (DG_FP)(p.second.values.size() - 1);
      s_variance.push_back(variance);
    }
    */
  }
}

/*
void PolyApprox::stencil_data_variance(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                              const DG_FP *y_ptr, const DG_FP *s_ptr, vector<DG_FP> &x,
                              vector<DG_FP> &y, vector<DG_FP> &s_mean, vector<DG_FP> &s_variance) {
  map<DGUtils::Vec<2>, Point> pointMap;
  get_stencil_values(cell_ind, stencil, x_ptr, y_ptr, s_ptr, offset_x, offset_y, pointMap);

  for(auto const &p : pointMap) {
    x.push_back(p.second.coord[0]);
    y.push_back(p.second.coord[1]);
    // Calculate mean
    DG_FP sum = 0.0;
    for(const auto &val : p.second.values) {
      sum += val;
    }
    const DG_FP mean = sum / (DG_FP)p.second.values.size();
    s_mean.push_back(mean);
    // Calculate variance
    if(p.second.values.size() == 1) {
      s_variance.push_back(0.0);
    } else {
      sum = 0.0;
      for(const auto &val : p.second.values) {
        sum += (val - mean) * (val - mean);
      }
      const DG_FP variance = sum / (DG_FP)(p.second.values.size() - 1);
      s_variance.push_back(variance);
    }
  }
}
*/

PolyApprox::PolyApprox(const int cell_ind, set<int> stencil,
                       const DG_FP *x_ptr, const DG_FP *y_ptr,
                       const DG_FP *s_ptr, const DG_FP h) {
  calc_offset(cell_ind, x_ptr, y_ptr);

  vector<DG_FP> x_vec, y_vec, s_mean_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, s_ptr, x_vec, y_vec, s_mean_vec);

  if(do_stencil_correction) {
    std::vector<int> bodies = body_scan(x_vec, y_vec, s_mean_vec);
    local_stencil_correction(x_vec, y_vec, s_mean_vec, bodies);
  }

  // Calc h
  auto min_x = std::min_element(x_vec.begin(), x_vec.end());
  auto min_y = std::min_element(y_vec.begin(), y_vec.end());
  auto max_x = std::max_element(x_vec.begin(), x_vec.end());
  auto max_y = std::max_element(y_vec.begin(), y_vec.end());

  DGUtils::Vec<2> min_pt(*min_x, *min_y);
  DGUtils::Vec<2> max_pt(*max_x, *max_y);

  DG_FP new_h = (max_pt - min_pt).magnitude();

  fit_poly(x_vec, y_vec, s_mean_vec, new_h);
}

PolyApprox::PolyApprox(std::vector<DG_FP> &c, DG_FP off_x, DG_FP off_y) {
  offset_x = off_x;
  offset_y = off_y;
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(c[i]);
  }
}

bool double_eq(const double d0, const double d1) {
  return fabs(d0 - d1) < 1e-8;
}

std::vector<std::set<int>> PolyApprox::construct_adj_list(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  std::vector<DGUtils::Vec<2>> pts;
  for(int i = 0; i < x.size(); i++) {
    pts.push_back(DGUtils::Vec<2>(x[i], y[i]));
  }

  // Delunay triangulation of DG nodes in stencil
  CDT::Triangulation<double> cdt;
  cdt.insertVertices(pts.begin(), pts.end(), [](const DGUtils::Vec<2>& p){ return p[0]; }, [](const DGUtils::Vec<2>& p){ return p[1]; });
  cdt.eraseSuperTriangle();
  auto triangles = cdt.triangles;
  auto vertices = cdt.vertices;

  std::map<int,int> vertices_map;
  for(int i = 0; i < pts.size(); i++) {
    for(int j = 0; j < pts.size(); j++) {
      if(double_eq(vertices[i].x, pts[j][0]) && double_eq(vertices[i].y, pts[j][1])) {
        vertices_map.insert({i, j});
        break;
      }
    }
  }

  // Construct adjaciency list from triangulation
  std::vector<std::set<int>> adj_map(pts.size());
  for(const auto &tri : triangles) {
    const int node_0 = vertices_map.at(tri.vertices[0]);
    const int node_1 = vertices_map.at(tri.vertices[1]);
    const int node_2 = vertices_map.at(tri.vertices[2]);

    adj_map[node_0].insert(node_1);
    adj_map[node_0].insert(node_2);
    adj_map[node_1].insert(node_0);
    adj_map[node_1].insert(node_2);
    adj_map[node_2].insert(node_0);
    adj_map[node_2].insert(node_1);
  }

  return adj_map;
}

std::vector<int> PolyApprox::body_scan(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s) {
  // Construct adjaciency list from triangulation
  std::vector<std::set<int>> adj_map = construct_adj_list(x, y);

  // Do the body scan
  std::vector<int> bodies(x.size(), -1);
  std::set<int> visited_pts;
  int body_ind = 0;
  int current_pt = 0;
  bool current_body_is_pos = s[current_pt] >= 0.0;
  bodies[current_pt] = body_ind;
  visited_pts.insert(current_pt);
  std::stack<int> query_queue;
  for(const auto &ind : adj_map[current_pt]) {
    if(s[ind] >= 0.0 == current_body_is_pos && visited_pts.count(ind) == 0) {
      bodies[ind] = body_ind;
      query_queue.push(ind);
    }
  }

  while(visited_pts.size() != x.size()) {
    if(query_queue.size() == 0) {
      body_ind++;
      for(int i = 0; i < x.size(); i++) {
        if(visited_pts.count(i) == 0) {
          current_pt = i;
          break;
        }
      }
      current_body_is_pos = s[current_pt] >= 0.0;
      bodies[current_pt] = body_ind;
      visited_pts.insert(current_pt);
      for(const auto &ind : adj_map[current_pt]) {
        if(s[ind] >= 0.0 == current_body_is_pos && visited_pts.count(ind) == 0) {
          bodies[ind] = body_ind;
          query_queue.push(ind);
        }
      }
    } else {
      current_pt = query_queue.top();
      query_queue.pop();
      visited_pts.insert(current_pt);
      for(const auto &ind : adj_map[current_pt]) {
        if(s[ind] >= 0.0 == current_body_is_pos && visited_pts.count(ind) == 0) {
          bodies[ind] = body_ind;
          query_queue.push(ind);
        }
      }
    }
  }

  return bodies;
}

DG_FP min_distance_to_line_segment(const DGUtils::Vec<2> &seg0, const DGUtils::Vec<2> &seg1, const DGUtils::Vec<2> &pt) {
  DG_FP mag_squared = (seg0 - seg1).sqr_magnitude();
  if(mag_squared < 1e-14)
    return (seg0 - pt).magnitude();

  DGUtils::Vec<2> tmp_0 = pt - seg0;
  DGUtils::Vec<2> tmp_1 = seg1 - seg0;
  DG_FP tmp = tmp_0.dot(tmp_1) / mag_squared;
  DG_FP factor = fmax(0, fmin(1.0, tmp));
  DGUtils::Vec<2> closest_pt = seg0 + factor * tmp_1;
  return (closest_pt - pt).magnitude();
}

bool would_make_cycle(const std::pair<int,int> &new_seg, const std::set<std::pair<int,int>> &segments, const int num_points) {
  std::vector<std::set<int>> adj_list(num_points);
  for(const auto &seg : segments) {
    adj_list[seg.first].insert(seg.second);
    adj_list[seg.second].insert(seg.first);
  }
  adj_list[new_seg.first].insert(new_seg.second);
  adj_list[new_seg.second].insert(new_seg.first);

  std::vector<bool> visited(num_points, false);
  int current_node = new_seg.first;
  int parent_node = -1;
  std::stack<int> backtrack;
  std::set<int> backtrack_contents;
  while(true) {
    visited[current_node] = true;
    int next_node = -1;
    for(const auto &ind : adj_list[current_node]) {
      if(visited[ind]) {
        if(ind != parent_node && backtrack_contents.count(ind) != 0)
          return true;
      } else {
        next_node = ind;
      }
    }

    // Handle backtracking
    if(next_node == -1) {
      if(backtrack.size() == 0) {
        return false;
      }

      current_node = backtrack.top();
      backtrack.pop();
      backtrack_contents.erase(current_node);
      if(backtrack.size() == 0) {
        parent_node = -1;
      } else {
        parent_node = backtrack.top();
      }
    } else {
      backtrack.push(current_node);
      backtrack_contents.insert(current_node);
      parent_node = current_node;
      current_node = next_node;
    }
  }
}

void PolyApprox::local_stencil_correction(std::vector<DG_FP> &x, std::vector<DG_FP> &y, std::vector<DG_FP> &s, std::vector<int> &bodies) {
  // TODO reuse adj_list from body_scan
  std::vector<std::set<int>> adj_list = construct_adj_list(x, y);

  // Count how many bodies there were
  int num_bodies = 0;
  for(int i = 0; i < x.size(); i++) {
    num_bodies = num_bodies < bodies[i] ? bodies[i] : num_bodies;
  }
  num_bodies++;

  if(num_bodies <= 2)
    return;

  if(num_bodies != 3) {
    printf("Number of bodies in stencil is %d\n", num_bodies);
    return;
  }

  // Find min distance from origin for each body
  DG_FP body_dists[3] = {1e5, 1e5, 1e5};
  std::set<int> pos_bodies, neg_bodies;
  for(int i = 0; i < x.size(); i++) {
    DG_FP dist = x[i] * x[i] + y[i] * y[i];
    body_dists[bodies[i]] = body_dists[bodies[i]] > dist ? dist : body_dists[bodies[i]];
    if(s[i] >= 0.0)
      pos_bodies.insert(bodies[i]);
    else
      neg_bodies.insert(bodies[i]);
  }

  // Decide which body is problematic
  int problem_body = -1;
  int safe_body = -1;
  if(pos_bodies.size() > neg_bodies.size()) {
    DG_FP max_dist = 0.0;
    for(const auto &body_ind : pos_bodies) {
      if(max_dist < body_dists[body_ind]) {
        max_dist = body_dists[body_ind];
        problem_body = body_ind;
      }
    }

    for(const auto &body_ind : pos_bodies) {
      if(body_ind != problem_body)
        safe_body = body_ind;
    }
  } else {
    DG_FP max_dist = 0.0;
    for(const auto &body_ind : neg_bodies) {
      if(max_dist < body_dists[body_ind]) {
        max_dist = body_dists[body_ind];
        problem_body = body_ind;
      }
    }

    for(const auto &body_ind : neg_bodies) {
      if(body_ind != problem_body)
        safe_body = body_ind;
    }
  }

  // Label nodes according to Ngo paper
  std::vector<int> node_labels(x.size(), -1);
  for(int i = 0; i < x.size(); i++) {
    if(bodies[i] == safe_body) {
      node_labels[i] = 1;
      for(const auto &ind : adj_list[i]) {
        if(node_labels[ind] == -1)
          node_labels[ind] = 2;
      }
    }
  }
  for(int i = 0; i < x.size(); i++) {
    if(node_labels[i] == -1)
      node_labels[i] = 3;
  }
/*
  // std::cout << "*** NODE LABELS ***" << std::endl;
  std::ofstream nodes_file("nodes.csv");
  nodes_file << "X,Y,S,B,NL" << std::endl;
  for(int i = 0; i < x.size(); i++) {
    nodes_file << x[i] << ",";
    nodes_file << y[i] << ",";
    nodes_file << s[i] << ",";
    nodes_file << bodies[i] << ",";
    nodes_file << node_labels[i] << std::endl;
  }
  nodes_file.close();
*/
  // Get points on interface using linear interpolation between nodes of label 1 and 2
  std::vector<DGUtils::Vec<2>> interface_pts;
  for(int i = 0; i < x.size(); i++) {
    if(node_labels[i] == 2) {
      for(const auto &ind : adj_list[i]) {
        if(node_labels[ind] == 1) {
          DGUtils::Vec<2> pt;
          DG_FP factor = s[i] / (s[i] - s[ind]);
          pt[0] = x[i] + factor * (x[ind] - x[i]);
          pt[1] = y[i] + factor * (y[ind] - y[i]);
          interface_pts.push_back(pt);
        }
      }
    }
  }
/*
  // std::cout << "*** INTERFACE PTS ***" << std::endl;
  std::ofstream pts_file("pts.csv");
  pts_file << "X,Y" << std::endl;
  for(int i = 0; i < interface_pts.size(); i++) {
    pts_file << interface_pts[i].x << ",";
    pts_file << interface_pts[i].y << std::endl;
  }
  pts_file.close();
*/
/*
  // Local reinit of nodes with label 3 just using points on interface
  bool node_3_neg = pos_bodies.size() > neg_bodies.size();
  for(int i = 0; i < x.size(); i++) {
    if(node_labels[i] == 3) {
      // Find closest point
      DG_FP dist = 1e5;
      for(int j = 0; j < interface_pts.size(); j++) {
        DG_FP tmp_dist = (interface_pts[j].x - x[i]) * (interface_pts[j].x - x[i])
                       + (interface_pts[j].y - y[i]) * (interface_pts[j].y - y[i]);
        dist = dist > tmp_dist ? tmp_dist : dist;
      }
      dist = sqrt(dist);
      s[i] = dist;
      if(node_3_neg) s[i] = -s[i];
    }
  }
*/

  // Get segments of line representing interface
  std::set<std::pair<int,int>> segments;
  const int total_number_of_segments = interface_pts.size() - 1;
  std::vector<int> pt_degree(interface_pts.size(), 0);
  while(segments.size() < total_number_of_segments) {
    // Iterate over all potential segments, find shortest one
    std::pair<int,int> shortest_segment = {-1, -1};
    DG_FP dist = 1e5;
    for(int i = 0; i < interface_pts.size(); i++) {
      if(pt_degree[i] == 2)
        continue;
      for(int j = i + 1; j < interface_pts.size(); j++) {
        if(pt_degree[j] == 2 || segments.count({i, j}) != 0 || would_make_cycle({i, j}, segments, interface_pts.size()))
          continue;
        DG_FP tmp_dist = (interface_pts[i] - interface_pts[j]).sqr_magnitude();
        if(tmp_dist < dist) {
          dist = tmp_dist;
          shortest_segment = {i, j};
        }
      }
    }

    pt_degree[shortest_segment.first]++;
    pt_degree[shortest_segment.second]++;
    segments.insert(shortest_segment);
  }
/*
  std::ofstream segments_file("segments.csv");
  segments_file << "X,Y" << std::endl;
  for(const auto &seg : segments) {
    for(int i = 0; i < 20; i++) {
      DG_FP factor = (DG_FP) i / 19.0;
      segments_file << interface_pts[seg.first].x + factor * (interface_pts[seg.second].x - interface_pts[seg.first].x) << ",";
      segments_file << interface_pts[seg.first].y + factor * (interface_pts[seg.second].y - interface_pts[seg.first].y) << std::endl;
    }
  }
  segments_file.close();
  exit(-1);
*/
  // Local reinit of nodes with label 3
  bool node_3_neg = pos_bodies.size() > neg_bodies.size();
  for(int i = 0; i < x.size(); i++) {
    if(node_labels[i] == 3) {
      DGUtils::Vec<2> current_pt(x[i], y[i]);
      // Find closest point
      DG_FP dist = 1e5;
      int closest_pt = -1;
      for(int j = 0; j < interface_pts.size(); j++) {
        DG_FP tmp_dist = (interface_pts[j] - current_pt).sqr_magnitude();
        if(dist > tmp_dist) {
          dist = tmp_dist;
          closest_pt = j;
        }
      }

      // Get line segments involving the closest point
      std::vector<std::pair<int,int>> segments_to_consider;
      for(const auto &seg : segments) {
        if(seg.first == closest_pt || seg.second == closest_pt)
          segments_to_consider.push_back(seg);
      }

      // Get smallest distance to the segements
      DG_FP reinit_dist = 1e5;
      for(const auto &seg : segments_to_consider) {
        DG_FP tmp = min_distance_to_line_segment(interface_pts[seg.first], interface_pts[seg.second], current_pt);
        dist = dist > tmp ? tmp : dist;
      }

      s[i] = dist;
      if(node_3_neg) s[i] = -s[i];
    }
  }

/*
  // std::cout << "*** REINIT NODES ***" << std::endl;
  std::ofstream reinit_nodes_file("reinit_nodes.csv");
  reinit_nodes_file << "X,Y,S,B,NL" << std::endl;
  for(int i = 0; i < x.size(); i++) {
    reinit_nodes_file << x[i] << ",";
    reinit_nodes_file << y[i] << ",";
    reinit_nodes_file << s[i] << ",";
    reinit_nodes_file << bodies[i] << ",";
    reinit_nodes_file << node_labels[i] << std::endl;
  }
  reinit_nodes_file.close();

  exit(-1);
*/
}

std::vector<int> PolyApprox::get_labels_for_least_squares(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s) {
  // Construct adjaciency list from triangulation
  std::vector<std::set<int>> adj_list = construct_adj_list(x, y);
  std::vector<int> bodies = body_scan(x, y, s);

  // Label nodes with distance to interface
  std::vector<int> labels(x.size(), -1);
  bool pts_not_visited = false;
  for(int i = 0; i < x.size(); i++) {
    const int current_body = bodies[i];
    for(const int &ind : adj_list[i]) {
      if(bodies[ind] != current_body)
        labels[i] = 0;
    }
    if(labels[i] == -1)
      pts_not_visited = true;
  }

  int current_label = 1;
  while(pts_not_visited) {
    pts_not_visited = false;
    for(int i = 0; i < x.size(); i++) {
      if(labels[i] != -1) continue;
      for(const int &ind : adj_list[i]) {
        if(labels[ind] == current_label - 1)
          labels[i] = current_label;
      }
      if(labels[i] == -1)
        pts_not_visited = true;
    }
    current_label++;
  }

  return labels;
}

arma::vec PolyApprox::get_weights_for_least_squares(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s,
                                                    const DG_FP h) {
  // std::vector<int> labels = get_labels_for_least_squares(x, y, s);
  // arma::vec w(x.size());
  // for(int i = 0; i < x.size(); i++) {
  //   w(i) = 1.0 / (1.0 + 10.0 * labels[i]);
  // }

  // arma::vec w(x.size());
  // for(int i = 0; i < x.size(); i++) {
  //   DG_FP sd = sqrt(fmax(1e-12, s_variance[i]));
  //   w(i) = 1.0 / (sd * 1e3);
  // }

  const DG_FP sigma = h * 0.01;
  arma::vec w(x.size());
  for(int i = 0; i < x.size(); i++) {
    const DG_FP dist_from_origin = x[i] * x[i] + y[i] * y[i];
    w(i) = exp(-dist_from_origin / sigma);
  }

  return w;
}

void PolyApprox::fit_poly(vector<DG_FP> &x, vector<DG_FP> &y, vector<DG_FP> &s, const DG_FP h) {
  // Set A vandermonde matrix and b
  arma::mat A;
  if(RDF) {
    rdf_eps = 512.0 / (h * h);
    A = get_rdf_vandermonde(x, y);
    arma::vec b(s);
    arma::vec ans = arma::inv(A) * b;
    coeff = arma::conv_to<vector<DG_FP>>::from(ans);
    return;
  } else {
    A = get_vandermonde(x, y);
  }
  arma::vec b(s);
  arma::vec w = get_weights_for_least_squares(x, y, s, h);
  // std::vector<int> labels = get_labels_for_least_squares(x, y, s);

  bool node_with_wrong_sign = false;
  set<int> problem_inds;
  int redo_counter = 0;
  const int max_redo = 5;
  arma::vec ans;
  do {
    for(const int &ind : problem_inds) {
      w(ind) = w(ind) * 2.0;
    }
    arma::mat W = arma::diagmat(w);
    // arma::mat Q, R;
    // arma::qr(Q, R, W*A);
    // ans = arma::solve(R, Q.t() * W * b);
    ans = arma::solve(W * A, W * b);

    arma::vec res = A * ans;
    node_with_wrong_sign = false;
    problem_inds.clear();
    for(int i = 0; i < x.size(); i++) {
      if(res(i) > 0.0 != b(i) > 0.0) {
        node_with_wrong_sign = true;
        problem_inds.insert(i);
      }
    }

    redo_counter++;
  } while(node_with_wrong_sign && redo_counter < max_redo);

  coeff = arma::conv_to<vector<DG_FP>>::from(ans);

  // if(redo_counter == max_redo) printf("Max redo\n");
  // if(node_with_wrong_sign) printf("Node with wrong sign\n");
}

DG_FP PolyApprox::gaussian(DG_FP r) {
  return exp(-(rdf_eps * r));
}

arma::mat PolyApprox::get_rdf_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  for(int i = 0; i < x.size(); i++) {
    rdf_stencil_pts.push_back(DGUtils::Vec<2>(x[i], y[i]));
  }

  // Set A vandermonde matrix and b
  arma::mat A(rdf_stencil_pts.size(), rdf_stencil_pts.size());
  for(int i = 0; i < rdf_stencil_pts.size(); i++) {
    for(int j = 0; j < rdf_stencil_pts.size(); j++) {
      DG_FP distance_to_i = (rdf_stencil_pts[j] - rdf_stencil_pts[i]).sqr_magnitude();
      A(i, j) = gaussian(distance_to_i);
    }
  }

  return A;
}

#define C_POLY_IND 0
#define X_POLY_IND 1
#define Y_POLY_IND 2
#define X2_POLY_IND 3
#define XY_POLY_IND 4
#define Y2_POLY_IND 5
#define X3_POLY_IND 6
#define X2Y_POLY_IND 7
#define Y2X_POLY_IND 8
#define Y3_POLY_IND 9
#define X4_POLY_IND 10
#define X3Y_POLY_IND 11
#define X2Y2_POLY_IND 12
#define Y3X_POLY_IND 13
#define Y4_POLY_IND 14

arma::mat PolyApprox::get_2nd_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  arma::mat A(x.size(), 6);
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)  = 1.0;
    A(i,X_POLY_IND)  = x[i];
    A(i,Y_POLY_IND)  = y[i];
    A(i,X2_POLY_IND) = x[i] * x[i];
    A(i,XY_POLY_IND) = x[i] * y[i];
    A(i,Y2_POLY_IND) = y[i] * y[i];
  }

  return A;
}

arma::mat PolyApprox::get_3rd_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  // Set A vandermonde matrix and b
  arma::mat A(x.size(), 10);
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)   = 1.0;
    A(i,X_POLY_IND)   = x[i];
    A(i,Y_POLY_IND)   = y[i];
    A(i,X2_POLY_IND)  = x[i] * x[i];
    A(i,XY_POLY_IND)  = x[i] * y[i];
    A(i,Y2_POLY_IND)  = y[i] * y[i];
    A(i,X3_POLY_IND)  = x[i] * x[i] * x[i];
    A(i,X2Y_POLY_IND) = x[i] * x[i] * y[i];
    A(i,Y2X_POLY_IND) = x[i] * y[i] * y[i];
    A(i,Y3_POLY_IND)  = y[i] * y[i] * y[i];
  }

  return A;
}

arma::mat PolyApprox::get_4th_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  arma::mat A(x.size(), 15);
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)    = 1.0;
    A(i,X_POLY_IND)    = x[i];
    A(i,Y_POLY_IND)    = y[i];
    A(i,X2_POLY_IND)   = x[i] * x[i];
    A(i,XY_POLY_IND)   = x[i] * y[i];
    A(i,Y2_POLY_IND)   = y[i] * y[i];
    A(i,X3_POLY_IND)   = x[i] * x[i] * x[i];
    A(i,X2Y_POLY_IND)  = x[i] * x[i] * y[i];
    A(i,Y2X_POLY_IND)  = x[i] * y[i] * y[i];
    A(i,Y3_POLY_IND)   = y[i] * y[i] * y[i];
    A(i,X4_POLY_IND)   = x[i] * x[i] * x[i] * x[i];
    A(i,X3Y_POLY_IND)  = x[i] * x[i] * x[i] * y[i];
    A(i,X2Y2_POLY_IND) = x[i] * x[i] * y[i] * y[i];
    A(i,Y3X_POLY_IND)  = x[i] * y[i] * y[i] * y[i];
    A(i,Y4_POLY_IND)   = y[i] * y[i] * y[i] * y[i];
  }

  return A;
}

arma::mat PolyApprox::get_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y) {
  arma::mat res;
  if(N == 2) {
    res = get_2nd_order_vandermonde(x, y);
  } else if(N == 3) {
    res = get_3rd_order_vandermonde(x, y);
  } else if(N == 4) {
    res = get_4th_order_vandermonde(x, y);
  }
  return res;
}

DG_FP PolyApprox::val_at_2nd(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  return res;
}

DG_FP PolyApprox::val_at_3rd(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[Y2X_POLY_IND] * x * y * y;
  res += coeff[Y3_POLY_IND] * y * y * y;
  return res;
}

DG_FP PolyApprox::val_at_4th(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[Y2X_POLY_IND] * x * y * y;
  res += coeff[Y3_POLY_IND] * y * y * y;
  res += coeff[X4_POLY_IND] * x * x * x * x;
  res += coeff[X3Y_POLY_IND] * x * x * x * y;
  res += coeff[X2Y2_POLY_IND] * x * x * y * y;
  res += coeff[Y3X_POLY_IND] * x * y * y * y;
  res += coeff[Y4_POLY_IND] * y * y * y * y;
  return res;
}

DG_FP PolyApprox::val_at_rdf(const DG_FP x, const DG_FP y) {
  DGUtils::Vec<2> pt(x, y);

  DG_FP result = 0.0;
  for(int i = 0; i < rdf_stencil_pts.size(); i++) {
    DG_FP dist = (pt - rdf_stencil_pts[i]).sqr_magnitude();
    result += coeff[i] * gaussian(dist);
  }

  return result;
}

DG_FP PolyApprox::val_at(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  if(RDF) {
    res = val_at_rdf(x, y);
  } else if(N == 2) {
    res = val_at_2nd(x, y);
  } else if(N == 3) {
    res = val_at_3rd(x, y);
  } else if(N == 4) {
    res = val_at_4th(x, y);
  }
  return res;
}

void PolyApprox::grad_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
}

void PolyApprox::grad_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += coeff[Y2X_POLY_IND] * y * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * x * y;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;
}

void PolyApprox::grad_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += coeff[Y2X_POLY_IND] * y * y;
  dx += 4.0 * coeff[X4_POLY_IND] * x * x * x;
  dx += 3.0 * coeff[X3Y_POLY_IND] * x * x * y;
  dx += 2.0 * coeff[X2Y2_POLY_IND] * x * y * y;
  dx += coeff[Y3X_POLY_IND] * y * y * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * x * y;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;
  dy += coeff[X3Y_POLY_IND] * x * x * x;
  dy += 2.0 * coeff[X2Y2_POLY_IND] * x * x * y;
  dy += 3.0 * coeff[Y3X_POLY_IND] * x * y * y;
  dy += 4.0 * coeff[Y4_POLY_IND] * y * y * y;
}

void PolyApprox::grad_at_rdf(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  DGUtils::Vec<2> pt(x, y);
  for(int i = 0; i < rdf_stencil_pts.size(); i++) {
    DG_FP dist = (pt - rdf_stencil_pts[i]).sqr_magnitude();
    DG_FP gauss = gaussian(dist);
    dx += -2.0 * rdf_eps * x * coeff[i] * gauss;
    dy += -2.0 * rdf_eps * y * coeff[i] * gauss;
  }
}

void PolyApprox::grad_at(const DG_FP x, const DG_FP y,
                         DG_FP &dx, DG_FP &dy) {
  if(RDF) {
    grad_at_rdf(x, y, dx, dy);
  } else if(N == 2) {
    grad_at_2nd(x, y, dx, dy);
  } else if(N == 3) {
    grad_at_3rd(x, y, dx, dy);
  } else if(N == 4) {
    grad_at_4th(x, y, dx, dy);
  }
}

void PolyApprox::hessian_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2 = 2.0 * coeff[X2_POLY_IND];
  dxy = coeff[XY_POLY_IND];
  dy2 = 2.0 * coeff[Y2_POLY_IND];
}

void PolyApprox::hessian_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2  = 2.0 * coeff[X2_POLY_IND];
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;

  dxy  = coeff[XY_POLY_IND];
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;

  dy2  = 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;
}

void PolyApprox::hessian_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2  = 2.0 * coeff[X2_POLY_IND];
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;
  dx2 += 12.0 * coeff[X4_POLY_IND] * x * x;
  dx2 += 6.0 * coeff[X3Y_POLY_IND] * x * y;
  dx2 += 2.0 * coeff[X2Y2_POLY_IND] * y * y;

  dxy  = coeff[XY_POLY_IND];
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;
  dxy += 3.0 * coeff[X3Y_POLY_IND] * x * x;
  dxy += 4.0 * coeff[X2Y2_POLY_IND] * x * y;
  dxy += 3.0 * coeff[Y3X_POLY_IND] * y * y;

  dy2  = 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;
  dy2 += 2.0 * coeff[X2Y2_POLY_IND] * x * x;
  dy2 += 6.0 * coeff[Y3X_POLY_IND] * x * y;
  dy2 += 12.0 * coeff[Y4_POLY_IND] * y * y;
}

void PolyApprox::hessian_at_rdf(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2 = 0.0;
  dxy = 0.0;
  dy2 = 0.0;

  DGUtils::Vec<2> pt(x, y);
  for(int i = 0; i < rdf_stencil_pts.size(); i++) {
    DG_FP dist = (pt - rdf_stencil_pts[i]).sqr_magnitude();
    DG_FP gauss = gaussian(dist);
    dx2 += 2.0 * coeff[i] * rdf_eps * (2.0 * rdf_eps * x * x - 1.0) * gauss;
    dxy += 4.0 * coeff[i] * rdf_eps * rdf_eps * x * y * gauss;
    dy2 += 2.0 * coeff[i] * rdf_eps * (2.0 * rdf_eps * y * y - 1.0) * gauss;
  }
}

void PolyApprox::hessian_at(const DG_FP x, const DG_FP y,
                            DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  if(RDF) {
    hessian_at_rdf(x, y, dx2, dxy, dy2);
  } else if(N == 2) {
    hessian_at_2nd(x, y, dx2, dxy, dy2);
  } else if(N == 3) {
    hessian_at_3rd(x, y, dx2, dxy, dy2);
  } else if(N == 4) {
    hessian_at_4th(x, y, dx2, dxy, dy2);
  }
}

int PolyApprox::num_coeff() {
  if(N == 2) {
    return 6;
  } else if(N == 3) {
    return 10;
  } else if(N == 4) {
    return 15;
  } else {
    return -1;
  }
}

int PolyApprox::num_pts() {
  if(N == 2) {
    return 26;
  } else if(N == 3) {
    return 24;
  } else if(N == 4) {
    return 25;
  } else {
    return 0;
  }
}

int PolyApprox::num_elem_stencil() {
  if(N == 2) {
    return 13;
  } else if(N == 3) {
    return 13;
  } else if(N == 4) {
    return 13;
  } else {
    return 0;
  }
}

DG_FP PolyApprox::get_coeff(int ind) {
  return coeff[ind];
}

void PolyApprox::get_offsets(DG_FP &x, DG_FP &y) {
  x = offset_x;
  y = offset_y;
}

struct stencil_query {
  int ind;
  set<int> central_inds;
};

map<int,set<int>> PolyApprox::get_stencils(const set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  return single_layer_stencils(central_inds, edge_map, x_ptr, y_ptr);
/*
  timer->startTimer("PolyApprox - get_stencils");
  map<int,set<int>> stencils;
  map<int,stencil_query> queryInds;
  const int num_elements = num_elem_stencil();

  for(const auto &ind : central_inds) {
    set<int> st;
    st.insert(ind);
    stencils.insert({ind, st});
    stencil_query sq;
    sq.ind = ind;
    sq.central_inds.insert(ind);
    queryInds.insert({ind, sq});
  }

  const int numEdges = edge_map->from->size;
  while(queryInds.size() > 0) {
    map<int,stencil_query> newQueryInds;

    // Iterate over each edge pair
    for(int i = 0; i < numEdges * 2; i++) {
      // Find if this cell ind is in the query inds
      auto it = queryInds.find(edge_map->map[i]);
      if(it != queryInds.end()) {
        if(i % 2 == 0) {
          // For each central ind associated with this query ind
          for(const auto &ind : it->second.central_inds) {
            auto stencil_it = stencils.find(ind);
            // Check if the other cell in this edge is already in the stencil for this central ind
            if(stencil_it->second.find(edge_map->map[i + 1]) == stencil_it->second.end()
               && stencil_it->second.size() < num_elements) {
              stencil_it->second.insert(edge_map->map[i + 1]);
              // If stencil is not full then add to next rounds query inds
              if(stencil_it->second.size() < num_elements) {
                stencil_query sq;
                sq.ind = edge_map->map[i + 1];
                auto res = newQueryInds.insert({edge_map->map[i + 1], sq});
                res.first->second.central_inds.insert(ind);
              }
            }
          }
        } else {
          // For each central ind associated with this query ind
          for(const auto &ind : it->second.central_inds) {
            auto stencil_it = stencils.find(ind);
            // Check if the other cell in this edge is already in the stencil for this central ind
            if(stencil_it->second.find(edge_map->map[i - 1]) == stencil_it->second.end()
               && stencil_it->second.size() < num_elements) {
              stencil_it->second.insert(edge_map->map[i - 1]);
              // If stencil is not full then add to next rounds query inds
              if(stencil_it->second.size() < num_elements) {
                stencil_query sq;
                sq.ind = edge_map->map[i - 1];
                auto res = newQueryInds.insert({edge_map->map[i - 1], sq});
                res.first->second.central_inds.insert(ind);
              }
            }
          }
        }
      }
    }

    queryInds = newQueryInds;
  }
  timer->endTimer("PolyApprox - get_stencils");

  return stencils;
*/
}

bool share_coords(const DG_FP *x_ptr, const DG_FP *y_ptr, const std::vector<DGUtils::Vec<2>> &nodes) {
  for(int i = 0; i < DG_NP; i++) {
    for(int n = 0; n < 3; n++) {
      bool xCmp = abs(x_ptr[i] - nodes[n][0]) < 1e-8;
      bool yCmp = abs(y_ptr[i] - nodes[n][1]) < 1e-8;
      if(xCmp && yCmp) return true;
    }
  }
  return false;
}

map<int,set<int>> PolyApprox::single_layer_stencils(const set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  timer->startTimer("PolyApprox - get_stencils");
  map<int,set<int>> stencils;
  map<int,stencil_query> queryInds;
  const int num_elements = num_elem_stencil();
  map<int,std::vector<DGUtils::Vec<2>>> central_inds_nodes;

  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];

  for(const auto &ind : central_inds) {
    set<int> st;
    st.insert(ind);
    stencils.insert({ind, st});
    stencil_query sq;
    sq.ind = ind;
    sq.central_inds.insert(ind);
    queryInds.insert({ind, sq});

    std::vector<DGUtils::Vec<2>> nodes;
    DGUtils::Vec<2> node0, node1, node2;
    node0[0] = x_ptr[ind * DG_NP + fmask_node_ind_0];
    node0[1] = y_ptr[ind * DG_NP + fmask_node_ind_0];
    nodes.push_back(node0);
    node1[0] = x_ptr[ind * DG_NP + fmask_node_ind_1];
    node1[1] = y_ptr[ind * DG_NP + fmask_node_ind_1];
    nodes.push_back(node1);
    node2[0] = x_ptr[ind * DG_NP + fmask_node_ind_2];
    node2[1] = y_ptr[ind * DG_NP + fmask_node_ind_2];
    nodes.push_back(node2);
    central_inds_nodes.insert({ind, nodes});
  }

  const int numEdges = edge_map->from->size;

  while(queryInds.size() > 0) {
    map<int,stencil_query> newQueryInds;

    // Iterate over each edge pair
    for(int i = 0; i < numEdges * 2; i++) {
      // Find if this cell ind is in the query inds
      auto it = queryInds.find(edge_map->map[i]);
      if(it != queryInds.end()) {
        if(i % 2 == 0) {
          // For each central ind associated with this query ind
          for(const auto &ind : it->second.central_inds) {
            auto stencil_it = stencils.find(ind);
            // Check if the other cell in this edge is already in the stencil for this central ind
            if(stencil_it->second.find(edge_map->map[i + 1]) == stencil_it->second.end()
               && stencil_it->second.size() < num_elements) {
              // Check if we share a node with the central ind
              auto node_coords = central_inds_nodes.at(ind);
              if(share_coords(x_ptr + edge_map->map[i + 1] * DG_NP, y_ptr + edge_map->map[i + 1] * DG_NP, node_coords)) {
                stencil_it->second.insert(edge_map->map[i + 1]);
                // If stencil is not full then add to next rounds query inds
                if(stencil_it->second.size() < num_elements) {
                  stencil_query sq;
                  sq.ind = edge_map->map[i + 1];
                  auto res = newQueryInds.insert({edge_map->map[i + 1], sq});
                  res.first->second.central_inds.insert(ind);
                }
              }
            }
          }
        } else {
          // For each central ind associated with this query ind
          for(const auto &ind : it->second.central_inds) {
            auto stencil_it = stencils.find(ind);
            // Check if the other cell in this edge is already in the stencil for this central ind
            if(stencil_it->second.find(edge_map->map[i - 1]) == stencil_it->second.end()
               && stencil_it->second.size() < num_elements) {
              auto node_coords = central_inds_nodes.at(ind);
              if(share_coords(x_ptr + edge_map->map[i - 1] * DG_NP, y_ptr + edge_map->map[i - 1] * DG_NP, node_coords)) {
                stencil_it->second.insert(edge_map->map[i - 1]);
                // If stencil is not full then add to next rounds query inds
                if(stencil_it->second.size() < num_elements) {
                  stencil_query sq;
                  sq.ind = edge_map->map[i - 1];
                  auto res = newQueryInds.insert({edge_map->map[i - 1], sq});
                  res.first->second.central_inds.insert(ind);
                }
              }
            }
          }
        }
      }
    }

    queryInds = newQueryInds;
  }
  timer->endTimer("PolyApprox - get_stencils");

  return stencils;
}

/*
 * Different ways of getting the stencil data
 */
void get_stencil_values(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                        const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP offset_x,
                        const DG_FP offset_y, map<DGUtils::Vec<2>, Point> &pointMap) {
  for(const auto &sten : stencil) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      DGUtils::Vec<2> coord;
      coord[0] = x_ptr[ind] - offset_x;
      coord[1] = y_ptr[ind] - offset_y;
      Point point;
      auto res = pointMap.insert(make_pair(coord, point));

      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.value = s_ptr[ind];
        res.first->second.count = 1;
      } else {
        // Point already exists
        res.first->second.value += s_ptr[ind];
        res.first->second.count++;
      }
    }
  }
}
