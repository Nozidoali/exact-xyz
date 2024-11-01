/*
 * Author: Hanyu Wang
 * Created time: 2024-03-31 14:45:41
 * Last Modified by: Hanyu Wang
 * Last Modified time: 2024-04-02 14:08:25
 */

#include "qgate.hpp"
#include "qstate.hpp"

#include <cmath>

namespace xyz
{
QState RY::operator()( const QState& state, const bool reverse ) const
{
  auto _theta = reverse ? -theta : theta;
  QState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    if ( new_state.index_to_weight.find( index ) == new_state.index_to_weight.end() )
      new_state.index_to_weight[index] = 0;
    new_state.index_to_weight[index] += std::cos( _theta / 2 ) * weight;
    uint32_t new_index = index ^ ( 1 << target );
    new_state.index_to_weight[new_index] += std::sin( _theta / 2 ) *
                                            ( ( index & ( 1 << target ) ) ? -weight : weight );
  }
  for ( auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end(); )
  {
    if ( std::abs( it->second ) < QState::eps )
      it = new_state.index_to_weight.erase( it );
    else
      ++it;
  }
  new_state.n_bits = state.n_bits;
  return new_state;
}
} // namespace xyz