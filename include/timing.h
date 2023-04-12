#ifndef __INS_TIMING_H
#define __INS_TIMING_H

#include <string>
#include <map>

class Timing {
public:
  void exportTimings(std::string filename, int iter, double time);

  void startTimer(const std::string &name);
  void endTimer(const std::string &name);
  double getTime(const std::string &name);

private:
  std::map<std::string,double> startTime, totalTime;
  std::map<std::string,int> numCalls;

  static const int col_width = 50;
};

#endif
