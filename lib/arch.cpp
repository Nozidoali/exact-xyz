#include "arch.hpp"

using json = nlohmann::json;
namespace xyz
{

namespace detail
{
Zone parseZone( const json& zone_json, Config& config )
{
  Zone zone;
  zone.zone_id = zone_json["zone_id"].get<int>();
  zone.offset = { zone_json["offset"][0].get<int>(), zone_json["offset"][1].get<int>() };
  zone.dimension = { zone_json["dimension"][0].get<int>(), zone_json["dimension"][1].get<int>() };

  for ( const auto& slm_json : zone_json["slms"] )
  {
    SLM slm;
    slm.id = slm_json["id"].get<int>();
    slm.site_separation = { slm_json["site_seperation"][0].get<int>(), slm_json["site_seperation"][1].get<int>() };
    slm.r = slm_json["r"].get<int>();
    slm.c = slm_json["c"].get<int>();
    slm.location = { slm_json["location"][0].get<int>(), slm_json["location"][1].get<int>() };
    zone.slm_ids.push_back( config.slms.size() );
    config.slm_id_to_index[slm.id] = config.slms.size();
    config.slms.push_back( slm );
  }
  return zone;
}
} // namespace detail

Config parseConfig( const std::string& filename )
{
  Config config;
  std::ifstream file( filename );
  json j;
  file >> j;

  // Parse name
  config.name = j["name"].get<std::string>();
  for ( const auto& zone_json : j["storage_zones"] )
    config.storage_zones.push_back( detail::parseZone( zone_json, config ) );
  for ( const auto& zone_json : j["entanglement_zones"] )
    config.entanglement_zones.push_back( detail::parseZone( zone_json, config ) );
  for ( const auto& aod_json : j["aods"] )
  {
    AOD aod;
    aod.id = aod_json["id"].get<int>();
    aod.site_separation = aod_json["site_seperation"].get<int>();
    aod.r = aod_json["r"].get<int>();
    aod.c = aod_json["c"].get<int>();
    config.aods.push_back( aod );
  }

  // Parse arch range
  config.arch_range = { { { j["arch_range"][0][0].get<int>(), j["arch_range"][0][1].get<int>() },
                          { j["arch_range"][1][0].get<int>(), j["arch_range"][1][1].get<int>() } } };

  // Parse rydberg range
  config.rydberg_range = { { { j["rydberg_range"][0][0][0].get<int>(), j["rydberg_range"][0][0][1].get<int>() },
                             { j["rydberg_range"][0][1][0].get<int>(), j["rydberg_range"][0][1][1].get<int>() } } };

  return config;
}

std::array<int, 2> Config::loc_slm( int slm_id, int r, int c ) const
{
  auto it = slm_id_to_index.find( slm_id );
  if ( it == slm_id_to_index.end() )
    throw std::runtime_error( "SLM ID not found" );
  auto slm = slms[it->second];
  return { slm.location[0] + c * slm.site_separation[0], slm.location[1] + r * slm.site_separation[1] };
}

} // namespace xyz