#include <cstdint>
#include <iostream>
#include <map>
#include <qcircuit.hpp>

using namespace xyz;
int main()
{
  QState state = dicke_state( 4, 2 );
  QCircuit qc = prepare_state( state );
  write_qasm2( qc, "temp.qasm" );
  std::cout << qc.num_cnots() << std::endl;
  QCircuit new_qc = resyn( qc );
  std::cout << new_qc.num_cnots() << std::endl;
  write_qasm2( new_qc, "temp_resyn.qasm" );
  return 0;
}