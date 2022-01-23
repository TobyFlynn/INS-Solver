#ifndef __INS_TIMING_H
#define __INS_TIMING_H

#include <string>

class Timing {
public:
  void exportTimings(std::string filename, int iter, double time);
  double getWallTime();
  double getMainLoop();

  void startWallTime();
  void endWallTime();

  void startSetup();
  void endSetup();

  void startMainLoop();
  void endMainLoop();

  void startAdvection();
  void endAdvection();

  void startPressure();
  void endPressure();

  void startViscosity();
  void endViscosity();

  void startSurface();
  void endSurface();

  void startSave();
  void endSave();

  void startPressureSetup();
  void endPressureSetup();

  void startPressureLinearSolve();
  void endPressureLinearSolve();

  void startViscositySetup();
  void endViscositySetup();

  void startViscosityLinearSolve();
  void endViscosityLinearSolve();

  void startKSPSolve();
  void endKSPSolve();

  void startBuildMat();
  void endBuildMat();

private:
  double cpu1, cpu2;

  double totalWallTime = 0.0;
  double wallTime1, wallTime2;

  double totalSetup = 0.0;
  double setup1, setup2;

  double totalMainLoop = 0.0;
  double mainLoop1, mainLoop2;

  double totalAdvection = 0.0;
  double advection1, advection2;

  double totalPressure = 0.0;
  double pressure1, pressure2;

  double totalViscosity = 0.0;
  double viscosity1, viscosity2;

  double totalSurface = 0.0;
  double surface1, surface2;

  double totalSave = 0.0;
  double save1, save2;

  double totalPressureSetup = 0.0;
  double pressureSetup1, pressureSetup2;

  double totalPressureLinearSolve = 0.0;
  double pressureLinearSolve1, pressureLinearSolve2;

  double totalViscositySetup = 0.0;
  double viscositySetup1, viscositySetup2;

  double totalViscosityLinearSolve = 0.0;
  double viscosityLinearSolve1, viscosityLinearSolve2;

  double totalKSPSolve = 0.0;
  double KSPSolve1, KSPSolve2;

  double totalBuildMat = 0.0;
  double buildMat1, buildMat2;
};

#endif
