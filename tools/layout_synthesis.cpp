#include <cmdline.hpp>
#include <cstdint>
#include <iostream>
#include <map>
#include <xyz.hpp>

using namespace xyz;
using namespace zoned;
using cmdline::parser;

parser CommandLineParser()
{
  parser opt;
  opt.add<std::string>( "arch", 'a', "path to the architecture JSON file", false, "../data/arch.json" );
  opt.add<std::string>( "input", 'i', "path to the input QASM2 file", false, "../data/input.qasm" );
  return opt;
}

int main( int argc, char** argv )
{
  auto opt = CommandLineParser();
  opt.parse_check( argc, argv );

  auto config = parseConfig( opt.get<std::string>( "arch" ) );
  auto qc = read_qasm2( opt.get<std::string>( "input" ) );

  // Output parsed data for verification
  std::cout << "Configuration Name: " << config.name << "\n";
  std::cout << "Number of Storage Zones: " << config.storage_zones.size() << "\n";
  std::cout << "Number of Entanglement Zones: " << config.entanglement_zones.size() << "\n";
  std::cout << "Number of AODs: " << config.aods.size() << "\n";

  // Synthesis
  layout_synthesis( qc, config );

  return 0;
}
