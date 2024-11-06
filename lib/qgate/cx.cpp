#include "qgate.hpp"

namespace xyz
{
QState CX::operator()( const QState& state, const bool reverse ) const
{
  (void)reverse; // the conjugate of CX is CX
  QState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    uint32_t new_index = index;
    if ( (bool)( ( index >> ctrl ) & (uint32_t)1u ) == phase )
      new_index ^= ( 1 << target );
    new_state.index_to_weight[new_index] = weight;
  }
  new_state.n_bits = state.n_bits;
  return new_state;
}
} // namespace xyz