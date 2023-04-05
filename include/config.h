#ifndef __INS_CONFIG_H
#define __INS_CONFIG_H

#include <string>

#include "inipp.h"

class Config {
public:
  Config(const std::string &filename);

  bool getStr(const std::string &section, const std::string &key, std::string &val);
  bool getInt(const std::string &section, const std::string &key, int &val);
  bool getDouble(const std::string &section, const std::string &key, double &val);

private:
  bool queryIni(const std::string &section, const std::string &key, std::string &val);

  inipp::Ini<char> ini;
};

#endif
