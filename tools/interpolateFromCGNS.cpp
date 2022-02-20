// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <getopt.h>
#include <cstring>
#include <algorithm>

#include "dg_utils.h"

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"points", required_argument, 0, 0},
  {"grid",   required_argument, 0, 0},
  {"sol",    required_argument, 0, 0},
  {0,    0,                  0,  0}
};

struct DGCell {
  double x[10];
  double y[10];
  double velx[10];
  double vely[10];
  double pr[10];
  int numPts;
};

struct cmpCoords {
    bool operator()(const pair<double,double>& a, const pair<double,double>& b) const {
        bool xCmp = abs(a.first - b.first) < 1e-8;
        bool yCmp = abs(a.second - b.second) < 1e-8;
        if(xCmp && yCmp) {
          return false;
        }
        if(xCmp) {
          return a.second < b.second;
        }
        return a.first < b.first;
    }
};

void load_points(const string &fileName, vector<double> &x, vector<double> &y);
void load_grid(const string &fileName, vector<double> &x, vector<double> &y, vector<int> &cells);
void load_sol(const string &fileName, vector<double> &x, vector<double> &y, vector<double> &velx, vector<double> &vely, vector<double> &pr);
void link_nodes_to_cells(const vector<double> &x, const vector<double> &y, const vector<int> &cell_map, vector<DGCell> &cells);
void link_sol_to_cells(const vector<double> &x, const vector<double> &y, const vector<double> &velx, const vector<double> &vely, const vector<double> &pr, vector<DGCell> &cells);
void interpolate_to_points(const vector<double> &x, const vector<double> &y, const vector<DGCell> &cells, vector<double> &velx, vector<double> &vely, vector<double> &pr);
void save_interpolated_points(const string &fileName, const vector<double> &x, const vector<double> &y, const vector<double> &velx, const vector<double> &vely, const vector<double> &pr);
bool is_point_in_cell(const double x, const double y, const DGCell &cell);

arma::vec r, s;

int main(int argc, char **argv) {
  string pointsFileName = "";
  string gridFileName = "";
  string solFileName = "";
  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"points") == 0) pointsFileName = optarg;
    if(strcmp((char*)options[opt_index].name,"grid") == 0) gridFileName = optarg;
    if(strcmp((char*)options[opt_index].name,"sol") == 0) solFileName = optarg;
  }

  arma::vec x_, y_;
  DGUtils::setRefXY(3, x_, y_);
  DGUtils::xy2rs(x_, y_, r, s);

  // Get points to interpolate to
  vector<double> x_points, y_points;
  load_points(pointsFileName, x_points, y_points);
  cout << "Loaded " << to_string(x_points.size()) << " from text file" << endl;

  // Read grid CGNS file
  vector<double> x_grid, y_grid;
  vector<int> cells_grid;
  load_grid(gridFileName, x_grid, y_grid, cells_grid);
  cout << "Loaded " << to_string(x_grid.size()) << " nodes and " << to_string(cells_grid.size() / 3) << " cells from grid CGNS file" << endl;

  // Convert to DGCell struct format
  vector<DGCell> cells;
  link_nodes_to_cells(x_grid, y_grid, cells_grid, cells);

  // Read in solution data from CGNS file
  vector<double> x_sol, y_sol, velx_sol, vely_sol, pr_sol;
  load_sol(solFileName, x_sol, y_sol, velx_sol, vely_sol, pr_sol);
  cout << "Loaded " << to_string(x_sol.size()) << " solution points" << endl;

  // Add solution data to DGCell struct
  link_sol_to_cells(x_sol, y_sol, velx_sol, vely_sol, pr_sol, cells);

  // Get solution data at new points
  vector<double> velx_points(x_points.size(), NAN);
  vector<double> vely_points(x_points.size(), NAN);
  vector<double> pr_points(x_points.size(), NAN);
  interpolate_to_points(x_points, y_points, cells, velx_points, vely_points, pr_points);

  // Save interpolated points to text file
  string outFile = "interpolated-points.txt";
  save_interpolated_points(outFile, x_points, y_points, velx_points, vely_points, pr_points);
  cout << "Saved interpolated points to \"" << outFile << "\"" << endl;

  return 0;
}

/*
 * Function to load points to interpolate to
 */
void load_points(const string &fileName, vector<double> &x, vector<double> &y) {
  fstream file;
  file.open(fileName, ios::in);
  if(!file.is_open()) {
    cerr << "ERROR: Could not open points text file" << endl << "Exiting" << endl;
    exit(-1);
  }

  // Read in points (each line contains one point "x y")
  string line;
  while(getline(file, line)) {
    int pos = line.find(" ");
    if(pos != string::npos) {
      string x_str = line.substr(0, pos);
      string y_str = line.substr(pos + 1, line.size());
      x.push_back(stod(x_str));
      y.push_back(stod(y_str));
    }
  }

  file.close();
}

/*
 * Helper functions for loading CGNS files
 */

template<typename T>
void cgns_load_cells(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells);

template<>
void cgns_load_cells<int>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  std::vector<cgsize_t> cells(numCells * 3);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  transform(cells.begin(), cells.end(), cgnsCells, [](int x) { return x - 1;});
}

template<>
void cgns_load_cells<long>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  std::vector<cgsize_t> cells(numCells * 3);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  transform(cells.begin(), cells.end(), cgnsCells, [](long x) { return (int)x - 1;});
}

/*
 * Function for loading grid CGNS file
 */
void load_grid(const string &fileName, vector<double> &x, vector<double> &y,
               vector<int> &cells) {
  // Read grid CGNS grid
  int file;
  if(cg_open(fileName.c_str(), CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, 1, 1, zoneName, &cg_numNodes);
  int numNodes = (int) cg_numNodes;

  // Get vertices
  x.resize(numNodes);
  y.resize(numNodes);
  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, 1, 1, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, 1, 1, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  // Get cell section
  char sectionName[33];
  CGNS_ENUMT(ElementType_t) elementType;
  cgsize_t elementStart, elementEnd;
  int elementNumBoundary, parentFlag;
  cg_section_read(file, 1, 1, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);

  // Get cells
  int numCells = elementEnd - elementStart + 1;
  cells.resize(numCells * 3);
  cgns_load_cells<cgsize_t>(file, 1, 1, cells.data(), numCells);

  cg_close(file);
}

/*
 * Function to load solution data from CGNS file
 */
void load_sol(const string &fileName, vector<double> &x, vector<double> &y,
              vector<double> &velx, vector<double> &vely, vector<double> &pr) {
  // Read grid CGNS grid
  int file;
  if(cg_open(fileName.c_str(), CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, 1, 1, zoneName, &cg_numNodes);
  int numNodes = (int) cg_numNodes;

  // Get vertices
  x.resize(numNodes);
  y.resize(numNodes);
  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, 1, 1, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, 1, 1, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  // Get solution data
  velx.resize(numNodes);
  vely.resize(numNodes);
  pr.resize(numNodes);
  cg_field_read(file, 1, 1, 1, "VelocityX", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, velx.data());
  cg_field_read(file, 1, 1, 1, "VelocityY", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, vely.data());
  cg_field_read(file, 1, 1, 1, "Pressure", CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, pr.data());

  cg_close(file);
}

/*
 * Function to generate the points of a cell
 */
void gen_pts(DGCell &cell) {
  const double n0x = cell.x[0];
  const double n0y = cell.y[0];
  const double n1x = cell.x[1];
  const double n1y = cell.y[1];
  const double n2x = cell.x[2];
  const double n2y = cell.y[2];

  arma::vec ones = arma::ones<arma::vec>(10);

  arma::vec x = 0.5 * ((-r - s) * n0x + (ones + r) * n1x + (ones + s) * n2x);
  arma::vec y = 0.5 * ((-r - s) * n0y + (ones + r) * n1y + (ones + s) * n2y);

  for(int i = 0; i < 10; i++) {
    cell.x[i] = x[i];
    cell.y[i] = y[i];
    cell.velx[i] = NAN;
    cell.vely[i] = NAN;
    cell.pr[i]   = NAN;
  }

  cell.numPts = 10;
}

/*
 * Function for linking nodes to a cell struct
 */
void link_nodes_to_cells(const vector<double> &x, const vector<double> &y,
                         const vector<int> &cell_map, vector<DGCell> &cells) {
  const int numCells = cell_map.size() / 3;
  for(int i = 0; i < numCells; i++) {
    DGCell c;
    for(int j = 0; j < 3; j++) {
      c.x[j] = x[cell_map[i * 3 + j]];
      c.y[j] = y[cell_map[i * 3 + j]];
    }
    c.numPts = 3;
    gen_pts(c);
    cells.push_back(c);
  }
}

/*
 * Function to link solution data to cells
 */
struct Point {
  double x;
  double y;
  vector<int> cells;
};

void add_sol_to_cell(const double x, const double y, const double velx,
                     const double vely, const double pr, DGCell &cell) {
  for(int i = 0; i < 10; i++) {
    if(abs(x - cell.x[i]) < 1e-8 && abs(y - cell.y[i]) < 1e-8) {
      cell.velx[i] = velx;
      cell.vely[i] = vely;
      cell.pr[i]   = pr;
      return;
    }
  }
  cerr << "ERROR: could not find point within cell to add solution to" << endl;
  exit(-1);
}

void link_sol_to_cells(const vector<double> &x, const vector<double> &y,
                       const vector<double> &velx, const vector<double> &vely,
                       const vector<double> &pr, vector<DGCell> &cells) {
  // Add all calculated grid points to a map and store pointers to cells that
  // each point is a part of
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < cells.size(); c++) {
    for(int i = 0; i < 10; i++) {
      pair<double,double> coords = make_pair(cells[c].x[i], cells[c].y[i]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      res.first->second->x = cells[c].x[i];
      res.first->second->y = cells[c].y[i];
      res.first->second->cells.push_back(c);
    }
  }

  // Now consider each solution point, get the cells that it is a part of and
  // add the solution to the DGCell struct
  for(int p = 0; p < x.size(); p++) {
    pair<double,double> coords = make_pair(x[p], y[p]);
    auto pt = pointMap.find(coords);
    if(pt == pointMap.end()) {
      cerr << "ERROR: Did not find solution point within point map" << endl;
      exit(-1);
    }
    int numCells = pt->second->cells.size();
    for(int i = 0; i < numCells; i++) {
      add_sol_to_cell(x[p], y[p], velx[p], vely[p], pr[p], cells[pt->second->cells[i]]);
    }
  }

  // Check for NAN values (i.e. does every point have solution data)
  int numNAN = 0;
  for(int i = 0; i < cells.size(); i++) {
    for(int j = 0; j < 10; j++) {
      if(isnan(cells[i].velx[j])) {
        numNAN++;
        cout << "(" << to_string(cells[i].x[j]) << "," << to_string(cells[i].y[j]) << ")" << endl;
      }
    }
  }

  if(numNAN == 0) {
    cout << "Solution values successfully added to all points" << endl;
  } else {
    cout << "ERROR: Some points do not have solution data associated with them. Number of NAN values: " << to_string(numNAN) << endl;
  }
}

/*
 * Function to interpolate solution to points
 */
void interpolate_to_points(const vector<double> &x, const vector<double> &y,
                           const vector<DGCell> &cells, vector<double> &velx,
                           vector<double> &vely, vector<double> &pr) {
  // Find cells that the points are within
  cout << "Searching through " << to_string(cells.size()) << " cells" << endl;
  for(int i = 0; i < cells.size(); i++) {
    for(int j = 0; j < x.size(); j++) {
      if(is_point_in_cell(x[j], y[j], cells[i])) {
        // cout << "Point (" << to_string(x[j]) << ", " << to_string(y[j]) << ") is in cell:" << endl;
        // cout << "  (" << to_string(cells[i].x[0]) << ", " << to_string(cells[i].y[0]) << ") (";
        // cout << to_string(cells[i].x[3]) << ", " << to_string(cells[i].y[3]) << ") (";
        // cout << to_string(cells[i].x[9]) << ", " << to_string(cells[i].y[9]) << ")" << endl;

        // Set new interp point at index 1
        arma::vec new_x(cells[i].x, 10);
        arma::vec new_y(cells[i].y, 10);
        new_x[1] = x[j];
        new_y[1] = y[j];

        // Calc new r and s points
        const double N0x = cells[i].x[0];
        const double N0y = cells[i].y[0];
        const double N1x = cells[i].x[3];
        const double N1y = cells[i].y[3];
        const double N2x = cells[i].x[9];
        const double N2y = cells[i].y[9];
        arma::vec new_s = (N1y - N0y) * (2.0 * new_x - N1x - N2x) - (N1x - N0x) * (2.0 * new_y - N1y - N2y);
        new_s = new_s / ((N0x - N1x) * (N2y - N0y) + (N1y - N0y) * (N2x - N0x));
        arma::vec new_r = 2.0 * new_x - new_s * (N2x - N0x) - N1x - N2x;
        new_r = new_r / (N1x - N0x);

        arma::mat V = DGUtils::vandermonde2D(r, s, 3);
        arma::mat invV = arma::inv(V);

        arma::mat newV = DGUtils::vandermonde2D(new_r, new_s, 3);
        arma::mat interpMat = newV * invV;

        arma::vec old_x(cells[i].x, 10);
        arma::vec old_y(cells[i].y, 10);
        arma::vec old_velx(cells[i].velx, 10);
        arma::vec old_vely(cells[i].vely, 10);
        arma::vec old_pr(cells[i].pr, 10);

        new_x = interpMat * old_x;
        new_y = interpMat * old_y;
        arma::vec new_velx = interpMat * old_velx;
        arma::vec new_vely = interpMat * old_vely;
        arma::vec new_pr = interpMat * old_pr;

        velx[j] = new_velx[1];
        vely[j] = new_vely[1];
        pr[j]   = new_pr[1];

        // cout << "New values for point (" << to_string(x[j]) << ", " << to_string(y[j]) << ") are:" << endl;
        // cout << "  (" << to_string(new_x[1]) << ", " << to_string(new_y[1]) << ") " << to_string(new_velx[1]) << "," << to_string(new_vely[1]) << "," << to_string(new_pr[1]) << endl;
      }
    }
  }
}

/*
 * Function to save interpolated points to a text file
 */
void save_interpolated_points(const string &fileName, const vector<double> &x,
                              const vector<double> &y, const vector<double> &velx,
                              const vector<double> &vely, const vector<double> &pr) {
  ofstream out(fileName.c_str());

  out << "X,Y,VelX,VelY,Pr" << endl;
  for(int i = 0; i < x.size(); i++) {
    out << to_string(x[i]) << "," << to_string(y[i]) << "," << to_string(velx[i]) << "," << to_string(vely[i]) << "," << to_string(pr[i]) << endl;
  }

  out.close();
}

/*
 * Function to check whether a point is contained with a triangular cell
 */
bool is_point_in_cell(const double x, const double y, const DGCell &cell) {
  double ABx = cell.x[3] - cell.x[0];
  double ABy = cell.y[3] - cell.y[0];

  double APx = x - cell.x[0];
  double APy = y - cell.y[0];

  double BCx = cell.x[9] - cell.x[3];
  double BCy = cell.y[9] - cell.y[3];

  double BPx = x - cell.x[3];
  double BPy = y - cell.y[3];

  double CAx = cell.x[0] - cell.x[9];
  double CAy = cell.y[0] - cell.y[9];

  double CPx = x - cell.x[9];
  double CPy = y - cell.y[9];

  double AB_AP = ABx * APy - ABy * APx;
  double BC_BP = BCx * BPy - BCy * BPx;
  double CA_CP = CAx * CPy - CAy * CPx;

  bool zero0 = AB_AP == 0.0 || AB_AP == -0.0;
  bool zero1 = BC_BP == 0.0 || BC_BP == -0.0;
  bool zero2 = CA_CP == 0.0 || CA_CP == -0.0;

  if(zero0 || zero1 || zero2) {
    // One zero means on an edge of the triangle
    // Two zeros means on a vertex of the triangle
    if(zero0 && !zero1 && !zero2) {
      return (BC_BP > 0.0 && CA_CP > 0.0) || (BC_BP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && zero1 && !zero2) {
      return (AB_AP > 0.0 && CA_CP > 0.0) || (AB_AP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && !zero1 && zero2) {
      return (AB_AP > 0.0 && BC_BP > 0.0) || (AB_AP < 0.0 && BC_BP < 0.0);
    } else if(zero0 && zero1 && !zero2) {
      return true;
    } else if(zero0 && !zero1 && zero2) {
      return true;
    } else if(!zero0 && zero1 && zero2) {
      return true;
    } else {
      cout << "ERROR: all zero case for point in triangle detection not implemented yet" << endl;
      exit(-1);
      return false;
    }
  }

  if(AB_AP > 0.0 && BC_BP > 0.0 && CA_CP > 0.0) {
    return true;
  } else if(AB_AP < 0.0 && BC_BP < 0.0 && CA_CP < 0.0) {
    return true;
  } else {
    return false;
  }
}
