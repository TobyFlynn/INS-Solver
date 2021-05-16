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

  void startLiftDrag();
  void endLiftDrag();

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

  void startLinearSolveMF();
  void endLinearSolveMF();

  void startLinearSolveMFMatMult();
  void endLinearSolveMFMatMult();

  void startLinearSolveMFRHS();
  void endLinearSolveMFRHS();

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

  double totalLiftDrag = 0.0;
  double liftDrag1, liftDrag2;

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

  double totalLinearSolveMF = 0.0;
  double linearSolveMF1, linearSolveMF2;

  double totalLinearSolveMFMatMult = 0.0;
  double linearSolveMFMatMult1, linearSolveMFMatMult2;

  double totalLinearSolveMFRHS = 0.0;
  double linearSolveMFRHS1, linearSolveMFRHS2;
};

#endif
