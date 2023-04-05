#include "config.h"

#include <fstream>

Config::Config(const std::string &filename) {
	std::ifstream is(filename);
	ini.parse(is);
  ini.strip_trailing_comments();
}

bool Config::getStr(const std::string &section, const std::string &key, std::string &val) {
  return queryIni(section, key, val);
}

bool Config::getInt(const std::string &section, const std::string &key, int &val) {
  std::string qry;
  if(queryIni(section, key, qry)) {
    val = stoi(qry);
    return true;
  } else {
    return false;
  }
}

bool Config::getDouble(const std::string &section, const std::string &key, double &val) {
  std::string qry;
  if(queryIni(section, key, qry)) {
    val = stod(qry);
    return true;
  } else {
    return false;
  }
}

bool Config::queryIni(const std::string &section, const std::string &key, std::string &val) {
  return inipp::get_value(ini.sections[section], key, val);
}
