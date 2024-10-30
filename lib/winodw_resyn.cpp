#include "qcircuit.hpp"

namespace xyz
{
namespace detail
{

template<bool verbose>
void resyn_impl( const QCircuit& circuit, QCircuit& new_circuit, const QState& initial_state, uint32_t pos )
{
  if ( pos >= circuit.pGates.size() )
    return;
  if constexpr ( verbose )
    std::cout << "initial_state: " << initial_state << std::endl;
  auto gate = circuit.pGates[pos];
  uint32_t target = gate->target;
  QState new_state = ( *gate )( initial_state );
  if constexpr ( verbose )
    std::cout << "gate: " << *gate << std::endl;
  uint32_t new_pos = pos + 1;
  for ( ; new_pos < circuit.pGates.size(); new_pos++ )
  {
    auto _gate = circuit.pGates[new_pos];
    if ( _gate->target != target )
      break;
    new_state = ( *_gate )( new_state );
    if constexpr ( verbose )
      std::cout << "gate: " << *_gate << std::endl;
  }
  if constexpr ( verbose )
    std::cout << "final_state: " << new_state << std::endl;
  resyn_impl<verbose>( circuit, new_circuit, new_state, new_pos );
}

} // namespace detail

QCircuit resyn( const QCircuit& circuit )
{
  QCircuit new_circuit( circuit.num_qbits );
  detail::resyn_impl<false>( circuit, new_circuit, ground_state( circuit.num_qbits ), 0 );
  return circuit;
}

} // namespace xyz