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
  opt.add<std::string>( "input", 'i', "path to the input JSON file", false, "../data/output.json" );
  return opt;
}

int main( int argc, char** argv )
{
  auto opt = CommandLineParser();
  opt.parse_check( argc, argv );

  auto config = parseConfig( opt.get<std::string>( "arch" ) );
  auto schedule = load_schedule( opt.get<std::string>( "input" ) );
  simulate( schedule, config );

  return 0;
}
