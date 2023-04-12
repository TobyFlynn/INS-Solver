#include "timing.h"

#include "op_seq.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <vector>

using namespace std;

#ifndef INS_MPI
void Timing::exportTimings(std::string filename, int iter, double time) {
  ofstream file(filename);

  file << std::left << std::setw(col_width) << "Iterations:" << iter << std::endl;
  file << std::left << std::setw(col_width) << "Final time:" << time << std::endl;
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << std::left << std::setw(col_width) << it->first + ":";
    file << std::left << std::setw(10) << it->second << numCalls[it->first] << std::endl;
  }

  file.close();
}
#else
#include "mpi.h"
void Timing::exportTimings(std::string filename, int iter, double time) {
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int num_timers = totalTime.size();

  // Check same number of timers on each process
  std::vector<int> num_timers_rcv(comm_size);
  MPI_Gather(&num_timers, 1, MPI_INT, num_timers_rcv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(rank == 0) {
    bool all_same = true;
    for(int i = 1; i < comm_size; i++) {
      if(num_timers_rcv[i] != num_timers)
        all_same = false;
    }

    if(!all_same) {
      // If not all the same then can't get all min max avg
      // So save timers on rank 0 then broadcast
      std::cout << "Could not do min max avg for timers, only saving timing data for rank 0" << std::endl;
      ofstream file(filename);
      file << std::left << std::setw(col_width) << "Iterations:" << iter << std::endl;
      file << std::left << std::setw(col_width) << "Final time:" << time << std::endl;
      for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
        file << std::left << std::setw(col_width) << it->first + ":" << it->second << std::endl;
      }
      file.close();

      int b = 0;
      MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return;
    } else {
      int b = 1;
      MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  } else {
    // Check if doing min max avg
    int b;
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(b == 0) {
      return;
    }
  }

  std::vector<double> times(num_timers);
  int i = 0;
  for(auto it = totalTime.begin(); it != totalTime.end(); it++) {
    times[i] = it->second;
    i++;
  }
  std::vector<double> avg(num_timers);
  std::vector<double> min(num_timers);
  std::vector<double> max(num_timers);

  MPI_Reduce(times.data(), avg.data(), num_timers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(times.data(), min.data(), num_timers, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(times.data(), max.data(), num_timers, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(rank != 0)
    return;

  // Get averages
  for(int j = 0; j < num_timers; j++) {
    avg[j] /= (double)comm_size;
  }

  ofstream file(filename);

  file << std::left << std::setw(col_width) << "Iterations:" << iter << std::endl;
  file << std::left << std::setw(col_width) << "Final time:" << time << std::endl;
  i = 0;
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << std::left << std::setw(col_width) << it->first + ":" << it->second;
    file << ",\t" << avg[i] << ",\t" << min[i] << ",\t" << max[i] << std::endl;
    i++;
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
  if(numCalls.count(name) == 0)
    numCalls[name] = 1;
  else
    numCalls[name]++;
}

double Timing::getTime(const std::string &name) {
  return totalTime.at(name);
}
