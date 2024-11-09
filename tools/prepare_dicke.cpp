#include <cstdint>
#include <fmt/format.h>
#include <iostream>
#include <map>
#include <xyz.hpp>

using namespace xyz;
int main()
{
  // QRState state_exp = dicke_state( 4, 2 );
  // QCircuit qc = prepare_state( state_exp );
  uint32_t min_qubits = 4;
  uint32_t max_qubits = 50;
  for ( uint32_t n_qbits = min_qubits; n_qbits < max_qubits; n_qbits++ )
  {
    for ( bool log_depth : { false, true } )
    {
      for ( bool opt_count : { false, true } )
      {
        fmt::print( "W state with {} qubits and log_depth = {}, opt_count = {}\n", n_qbits, log_depth, opt_count );
        QCircuit qc = prepare_w( n_qbits, log_depth, opt_count );
        fmt::print( "Number of CNOTs: {}\n", qc.num_cnots() );
        fmt::print( "Number of levels: {}\n", qc.lev_cnots() );
        qc = decompose_circuit( qc );
        int setup_idx = ( log_depth ? 1 : 0 ) + ( opt_count ? 2 : 0 );
        write_qasm2( qc, fmt::format( "w_{}_{}.qasm", n_qbits, setup_idx ) );
      }
    }
  }
  return 0;
}