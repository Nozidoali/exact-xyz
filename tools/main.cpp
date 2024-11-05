#include <arch.hpp>
#include <cstdint>
#include <iostream>
#include <map>
#include <qcircuit.hpp>

using namespace xyz;
int main()
{

  QState state = dicke_state( 5, 2 );
  QCircuit qc = prepare_state( state );
  qc = decompose_circuit( qc );
  write_qasm2( qc, "temp_decomposed.qasm" );

  QCircuit new_qc = resyn( qc );
  std::cout << new_qc.num_cnots() << std::endl;
  new_qc = decompose_circuit( new_qc );
  write_qasm2( new_qc, "temp_resyn.qasm" );
  return 0;
}