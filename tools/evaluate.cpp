#include <cmdline.hpp>
#include <cstdint>
#include <iostream>
#include <map>
#include <synthesis.hpp>

using namespace xyz;
using namespace xyz;
using cmdline::parser;

parser CommandLineParser()
{
  parser opt;
  opt.add<std::string>( "input", 'i', "path to the input JSON file", false, "../data/output.json" );
  return opt;
}

int main( int argc, char** argv )
{
  auto opt = CommandLineParser();
  opt.parse_check( argc, argv );

  auto schedule = load_schedule( opt.get<std::string>( "input" ) );
  simulate( schedule );

  return 0;
}
