#include "json.hpp"
#include "synthesis.hpp"

using json = nlohmann::json;

namespace xyz
{

// Main function to parse JSON
Schedule load_schedule( const std::string& filename )
{
  std::ifstream json_file( filename );
  json json_object;
  json_file >> json_object;

  std::string name = json_object["name"];
  std::string architecture_spec_path = json_object["architecture_spec_path"];
  Schedule schedule;

  for ( const auto& inst : json_object["instructions"] )
  {
    std::string type = inst["type"];
    if ( type == "init" )
    {
      int id = inst["id"];
      double begin_time = inst["begin_time"];
      double end_time = inst["end_time"];
      std::vector<std::vector<int>> init_locs;
      for ( const auto& loc : inst["init_locs"] )
      {
        init_locs.push_back( loc.get<std::vector<int>>() );
      }
      InitInst init_inst( id, begin_time, end_time, init_locs );
      schedule.push_back( std::make_shared<InitInst>( init_inst ) );
    }
    else if ( type == "1qGate" )
    {
      int id = inst["id"];
      std::string unitary = inst["unitary"];
      double begin_time = inst["begin_time"];
      double end_time = inst["end_time"];
      std::vector<std::vector<int>> locs;
      for ( const auto& loc : inst["locs"] )
      {
        locs.push_back( loc.get<std::vector<int>>() );
      }
      std::vector<std::pair<std::string, int>> gates;
      for ( const auto& gate : inst["gates"] )
      {
        gates.push_back( { gate["name"], gate["q"] } );
      }
      std::map<std::string, std::vector<int>> dependency;

      OneQGateInst oneq_inst( id, unitary, begin_time, end_time, locs, gates, dependency );
      schedule.push_back( std::make_shared<OneQGateInst>( oneq_inst ) );
    }
    else if ( type == "rearrangeJob" )
    {
      int id = inst["id"];
      int aod_id = inst["aod_id"];
      double begin_time = inst["begin_time"];
      double end_time = inst["end_time"];
      std::vector<std::vector<int>> begin_locs;
      for ( const auto& loc : inst["begin_locs"] )
      {
        begin_locs.push_back( loc.get<std::vector<int>>() );
      }
      std::vector<std::vector<int>> end_locs;
      for ( const auto& loc : inst["end_locs"] )
      {
        end_locs.push_back( loc.get<std::vector<int>>() );
      }
      std::vector<int> aod_qubits = inst["aod_qubits"].get<std::vector<int>>();
      std::map<std::string, std::vector<int>> dependency;

      RearrangeJobInst rearrange_inst( id, aod_id, begin_time, end_time, aod_qubits, begin_locs, end_locs, dependency );
      schedule.push_back( std::make_shared<RearrangeJobInst>( rearrange_inst ) );
      std::cout << "RearrangeJob instruction added " << rearrange_inst.id << std::endl;
    }
  }

  return schedule;
}

void simulate( const Schedule& schedule )
{
  for ( const auto& inst : schedule )
  {
    std::cout << "Printing instruction: ";
    inst->print();
  }
}

} // namespace xyz