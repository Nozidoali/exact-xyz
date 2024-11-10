#pragma once

#include "json.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace xyz {

namespace constants {
constexpr double t_transfer = 15;   // us
constexpr double t_rydberg  = 0.36; // us
constexpr double t_gate1q   = 52.0; // us
} // namespace constants

struct SLM {
  int id;
  std::array<int, 2> site_separation; // [row separation, column separation]
  int r;                              // Number of rows
  int c;                              // Number of columns
  std::array<int, 2> location;        // [x, y]
};

struct Zone {
  int zone_id;
  std::vector<uint32_t> slm_ids; // List of SLMs
  std::array<int, 2> offset;     // [x, y]
  std::array<int, 2> dimension;  // [width, height]
};

struct AOD {
  int id;
  int site_separation; // Site separation
  int r;               // Number of rows
  int c;               // Number of columns
};

class Config {
public:
  int config_id; // Configuration ID
  std::string name;
  std::vector<Zone> storage_zones;
  std::vector<Zone> entanglement_zones;
  std::vector<AOD> aods;
  std::vector<SLM> slms;
  std::array<std::array<int, 2>, 2> arch_range;    // [[x_min, y_min], [x_max, y_max]]
  std::array<std::array<int, 2>, 2> rydberg_range; // [[x_min, y_min], [x_max, y_max]]
  std::unordered_map<int, uint32_t> slm_id_to_index;

public:
  Config() = default;
  Config( int config_id, const std::string& name ) : config_id( config_id ), name( name ){};
  std::array<int, 2> loc_slm( int slm_id, int r, int c ) const;
};

Config parseConfig( const std::string& filename );

} // namespace xyz