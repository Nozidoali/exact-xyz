#pragma once

#include "json.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace xyz
{

namespace constants
{
  constexpr double t_transfer = 15; // us
  constexpr double t_rydberg = 0.36; // us
  constexpr double t_gate1q = 52.0; // us
} // namespace constants

struct SLM
{
  int id;
  std::array<int, 2> site_separation; // [row separation, column separation]
  int r;                              // Number of rows
  int c;                              // Number of columns
  std::array<int, 2> location;        // [x, y]
};

struct Zone
{
  int zone_id;
  std::vector<SLM> slms;        // List of SLMs
  std::array<int, 2> offset;    // [x, y]
  std::array<int, 2> dimension; // [width, height]
};

struct AOD
{
  int id;
  int site_separation; // Site separation
  int r;               // Number of rows
  int c;               // Number of columns
};

struct Config
{
  int config_id; // Configuration ID
  std::string name;
  std::vector<Zone> storage_zones;
  std::vector<Zone> entanglement_zones;
  std::vector<AOD> aods;
  std::array<std::array<int, 2>, 2> arch_range;    // [[x_min, y_min], [x_max, y_max]]
  std::array<std::array<int, 2>, 2> rydberg_range; // [[x_min, y_min], [x_max, y_max]]
};

Config parseConfig( const std::string& filename );

} // namespace xyz