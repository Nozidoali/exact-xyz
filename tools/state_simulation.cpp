#include <cmdline.hpp>
#include <cstdint>
#include <iostream>
#include <qcircuit.hpp>

using namespace xyz;
using cmdline::parser;

parser CommandLineParser()
{
  parser opt;
  opt.add<std::string>( "input", 'i', "path to the input QASM2 file", false, "../data/input.qasm" );
  return opt;
}

int main( int argc, char** argv )
{
  auto opt = CommandLineParser();
  opt.parse_check( argc, argv );

  auto qc = read_qasm2( opt.get<std::string>( "input" ) );

  uint32_t n_qbits = qc.num_qbits;
  QRState state = ground_rstate( n_qbits );
  QRState new_state = simulate_circuit( qc, state, true );

  std::cout << "Initial State: " << state << std::endl;
  std::cout << "Final State: " << new_state << std::endl;

  return 0;
}
