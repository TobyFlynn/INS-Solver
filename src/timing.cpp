#include "timing.h"

#include "op_seq.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

using namespace std;

#ifndef INS_MPI
void Timing::exportTimings(std::string filename, int iter, double time) {
  ofstream file(filename);

  file << std::left << std::setw(col_width) << "Iterations:" << iter << std::endl;
  file << std::left << std::setw(col_width) << "Final time:" << time << std::endl;
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << std::left << std::setw(col_width) << it->first + ":" << it->second << std::endl;
  }

  file.close();
}
#else
#include "mpi.h"
void Timing::exportTimings(std::string filename, int iter, double time) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank != 0)
    return;

  ofstream file(filename);

  file << std::left << std::setw(col_width) << "Iterations:" << iter << std::endl;
  file << std::left << std::setw(col_width) << "Final time:" << time << std::endl;
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << std::left << std::setw(col_width) << it->first + ":" << it->second << std::endl;
  }

  file.close();
}
#endif

void Timing::startTimer(const std::string &name) {
  double wall;
  struct timeval t;
  gettimeofday(&t, (struct timezone *)0);
  wall = t.tv_sec + t.tv_usec * 1.0e-6;
  startTime[name] = wall;
}

void Timing::endTimer(const std::string &name) {
  double wall;
  struct timeval t;
  gettimeofday(&t, (struct timezone *)0);
  wall = t.tv_sec + t.tv_usec * 1.0e-6;
  totalTime[name] += wall - startTime[name];
}

double Timing::getTime(const std::string &name) {
  return totalTime.at(name);
}
