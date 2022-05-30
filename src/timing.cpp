#include "timing.h"

#include "op_seq.h"

#include <fstream>
#include <iostream>

using namespace std;

#ifndef INS_MPI
void Timing::exportTimings(std::string filename, int iter, double time) {
  ofstream file(filename);

  // Create first row of csv file (headers of columns)
  file << "Iterations" << ",";
  file << "Final Time" << ",";
  file << "Wall Time" << ",";
  file << "Setup" << ",";
  file << "Main Loop" << ",";
  file << "Advection" << ",";
  file << "Pressure" << ",";
  file << "Viscosity" << ",";
  file << "Surface" << ",";
  file << "Save" << ",";
  file << "Pressure Setup" << ",";
  file << "Pressure Linear Solve" << ",";
  file << "Viscosity Setup" << ",";
  file << "Viscosity Linear Solve" << ",";
  file << "KSP Solve" << ",";
  file << "Build Mat" << ",";
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << it->first << ",";
  }
  file << endl;

  // Output timing data
  file << to_string(iter) << ",";
  file << to_string(time) << ",";
  file << to_string(totalWallTime) << ",";
  file << to_string(totalSetup) << ",";
  file << to_string(totalMainLoop) << ",";
  file << to_string(totalAdvection) << ",";
  file << to_string(totalPressure) << ",";
  file << to_string(totalViscosity) << ",";
  file << to_string(totalSurface) << ",";
  file << to_string(totalSave) << ",";
  file << to_string(totalPressureSetup) << ",";
  file << to_string(totalPressureLinearSolve) << ",";
  file << to_string(totalViscositySetup) << ",";
  file << to_string(totalViscosityLinearSolve) << ",";
  file << to_string(totalKSPSolve) << ",";
  file << to_string(totalBuildMat) << ",";
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << it->second << ",";
  }
  file << endl;

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

  // Create first row of csv file (headers of columns)
  file << "Iterations" << ",";
  file << "Final Time" << ",";
  file << "Wall Time" << ",";
  file << "Setup" << ",";
  file << "Main Loop" << ",";
  file << "Advection" << ",";
  file << "Pressure" << ",";
  file << "Viscosity" << ",";
  file << "Surface" << ",";
  file << "Save" << ",";
  file << "Pressure Setup" << ",";
  file << "Pressure Linear Solve" << ",";
  file << "Viscosity Setup" << ",";
  file << "Viscosity Linear Solve" << ",";
  file << "KSP Solve" << ",";
  file << "Build Mat" << ",";
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << it->first << ",";
  }
  file << endl;

  // Output timing data
  file << to_string(iter) << ",";
  file << to_string(time) << ",";
  file << to_string(totalWallTime) << ",";
  file << to_string(totalSetup) << ",";
  file << to_string(totalMainLoop) << ",";
  file << to_string(totalAdvection) << ",";
  file << to_string(totalPressure) << ",";
  file << to_string(totalViscosity) << ",";
  file << to_string(totalSurface) << ",";
  file << to_string(totalSave) << ",";
  file << to_string(totalPressureSetup) << ",";
  file << to_string(totalPressureLinearSolve) << ",";
  file << to_string(totalViscositySetup) << ",";
  file << to_string(totalViscosityLinearSolve) << ",";
  file << to_string(totalKSPSolve) << ",";
  file << to_string(totalBuildMat) << ",";
  for (auto it = totalTime.begin(); it != totalTime.end(); it++) {
    file << it->second << ",";
  }
  file << endl;

  file.close();
}
#endif

void Timing::startTimer(std::string name) {
  double cpu, wall;
  op_timers(&cpu, &wall);
  startTime[name] = wall;
}

void Timing::endTimer(std::string name) {
  double cpu, wall;
  op_timers(&cpu, &wall);
  totalTime[name] += wall - startTime[name];
}

double Timing::getWallTime() {
  return totalWallTime;
}

double Timing::getMainLoop() {
  return totalMainLoop;
}

void Timing::startWallTime() {
  op_timers(&cpu1, &wallTime1);
}
void Timing::endWallTime() {
  op_timers(&cpu2, &wallTime2);
  totalWallTime += wallTime2 - wallTime1;
}

void Timing::startSetup() {
  op_timers(&cpu1, &setup1);
}
void Timing::endSetup() {
  op_timers(&cpu2, &setup2);
  totalSetup += setup2 - setup1;
}

void Timing::startMainLoop() {
  op_timers(&cpu1, &mainLoop1);
}
void Timing::endMainLoop() {
  op_timers(&cpu2, &mainLoop2);
  totalMainLoop += mainLoop2 - mainLoop1;
}

void Timing::startAdvection() {
  op_timers(&cpu1, &advection1);
}
void Timing::endAdvection() {
  op_timers(&cpu2, &advection2);
  totalAdvection += advection2 - advection1;
}

void Timing::startPressure() {
  op_timers(&cpu1, &pressure1);
}
void Timing::endPressure() {
  op_timers(&cpu2, &pressure2);
  totalPressure += pressure2 - pressure1;
}

void Timing::startViscosity() {
  op_timers(&cpu1, &viscosity1);
}
void Timing::endViscosity() {
  op_timers(&cpu2, &viscosity2);
  totalViscosity += viscosity2 - viscosity1;
}

void Timing::startSurface() {
  op_timers(&cpu1, &surface1);
}
void Timing::endSurface() {
  op_timers(&cpu2, &surface2);
  totalSurface += surface2 - surface1;
}

void Timing::startSave() {
  op_timers(&cpu1, &save1);
}
void Timing::endSave() {
  op_timers(&cpu2, &save2);
  totalSave += save2 - save1;
}

void Timing::startPressureSetup() {
  op_timers(&cpu1, &pressureSetup1);
}
void Timing::endPressureSetup() {
  op_timers(&cpu2, &pressureSetup2);
  totalPressureSetup += pressureSetup2 - pressureSetup1;
}

void Timing::startPressureLinearSolve() {
  op_timers(&cpu1, &pressureLinearSolve1);
}
void Timing::endPressureLinearSolve() {
  op_timers(&cpu2, &pressureLinearSolve2);
  totalPressureLinearSolve += pressureLinearSolve2 - pressureLinearSolve1;
}

void Timing::startViscositySetup() {
  op_timers(&cpu1, &viscositySetup1);
}
void Timing::endViscositySetup() {
  op_timers(&cpu2, &viscositySetup2);
  totalViscositySetup += viscositySetup2 - viscositySetup1;
}

void Timing::startViscosityLinearSolve() {
  op_timers(&cpu1, &viscosityLinearSolve1);
}
void Timing::endViscosityLinearSolve() {
  op_timers(&cpu2, &viscosityLinearSolve2);
  totalViscosityLinearSolve += viscosityLinearSolve2 - viscosityLinearSolve1;
}

void Timing::startKSPSolve() {
  op_timers(&cpu1, &KSPSolve1);
}
void Timing::endKSPSolve() {
  op_timers(&cpu2, &KSPSolve2);
  totalKSPSolve += KSPSolve2 - KSPSolve1;
}

void Timing::startBuildMat() {
  op_timers(&cpu1, &buildMat1);
}
void Timing::endBuildMat() {
  op_timers(&cpu2, &buildMat2);
  totalBuildMat += buildMat2 - buildMat1;
}
