#include "measurements/2d/measurement_2d.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

Measurement2D::Measurement2D(INSSolverBase2D *i, const int sample_iter) {
  ins = i;
  sample_rate = sample_iter;
  count = 0;
  io_count = 0;
}

std::string Measurement2D::double_to_text(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}

bool Measurement2D::sample_this_iter() {
  count++;
  if(count != sample_rate) {
    return false;
  }
  count = 0;
  return true;
}

void Measurement2D::reset_io() {
  io_count = 0;
}

#ifndef DG_MPI
void Measurement2D::output(std::string &path) {
  std::ofstream file(path + get_filename() + ".txt");

  file << get_csv_header() << std::endl;
  std::string next_line = get_next_csv_line();
  while(next_line != "") {
    file << next_line << std::endl;
    next_line = get_next_csv_line();
  }

  file.close();
  reset_io();
}
#else
#include "mpi.h"
void Measurement2D::output(std::string &path) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
    std::ofstream file(path + get_filename() + ".txt");

    file << get_csv_header() << std::endl;
    std::string next_line = get_next_csv_line();
    while(next_line != "") {
      file << next_line << std::endl;
      next_line = get_next_csv_line();
    }

    file.close();
  }
  reset_io();
}
#endif