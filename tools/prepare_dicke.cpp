#include <cstdint>
#include <iostream>
#include <map>
#include <xyz.hpp>

using namespace xyz;
int main()
{
  QRState state_exp = dicke_state( 4, 2 );
  QCircuit qc = prepare_state( state_exp );
  qc = resyn( qc );
  // qc = decompose_circuit( qc );
  QRState state_act = simulate_circuit( qc, ground_rstate( 4 ), true );
  std::cout << state_act << std::endl;
  std::cout << state_exp << std::endl;
  write_qasm2( qc, "temp.qasm" );
  return 0;
}