#include <cstdint>
#include <fmt/format.h>
#include <iostream>
#include <map>
#include <cmdline.hpp>
#include <xyz.hpp>

using namespace xyz;
using namespace cmdline;

parser CommandLineParser()
{
  parser opt;
  opt.add<int>( "n", 'n', "number of qubits", true );
  opt.add<int>( "k", 'k', "number of excitations", true );
  return opt;
}

int main( int argc, char** argv )
{
  auto opt = CommandLineParser();
  opt.parse_check( argc, argv );
  auto n = opt.get<int>( "n" );
  auto k = opt.get<int>( "k" );
  auto circuit = prepare_dicke_state( n, k );
  write_qasm2( circuit, fmt::format( "dicke_{}_{}.qasm", n, k ) );
  (void)simulate_circuit( circuit, ground_rstate( n ), true );
  return 0;
}
