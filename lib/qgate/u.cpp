#include "qgate.hpp"
#include "qstate.hpp"

#include <cmath>

namespace xyz
{
bool Rotation::is_trivial( double theta, bool use_x )
{
  bool is_zero = std::abs( theta ) < eps || std::abs( theta - 2 * M_PI ) < eps;
  bool is_pi = std::abs( theta - M_PI ) < eps || std::abs( theta + M_PI ) < eps;
  if ( use_x )
    return is_zero || is_pi;
  return is_zero;
}

QRState X::operator()( const QRState& state, const bool reverse ) const
{
  (void)reverse;
  QRState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    uint32_t new_index = index ^ ( 1 << target );
    new_state.index_to_weight[new_index] = weight;
  }
  new_state.n_bits = state.n_bits;
  return new_state;
}

QRState RU2::operator()( const QRState& state, uint32_t target, const bool reverse ) const
{
  QRState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    if ( new_state.index_to_weight.find( index ) == new_state.index_to_weight.end() )
      new_state.index_to_weight[index] = 0;
    uint32_t new_index = index ^ ( 1 << target );
    if ( index & ( 1 << target ) )
    {
      new_state.index_to_weight[new_index] += c10[reverse] * weight;
      new_state.index_to_weight[index] += c11[reverse] * weight;
    }
    else
    {
      new_state.index_to_weight[index] += c00[reverse] * weight;
      new_state.index_to_weight[new_index] += c01[reverse] * weight;
    }
  }
  for ( auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end(); )
  {
    if ( std::abs( it->second ) < QRState::eps )
      it = new_state.index_to_weight.erase( it );
    else
      ++it;
  }
  new_state.n_bits = state.n_bits;
  return new_state;
}

QState U2::operator()( const QState& state, uint32_t target, const bool reverse ) const
{
  QState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    if ( new_state.index_to_weight.find( index ) == new_state.index_to_weight.end() )
      new_state.index_to_weight[index] = 0;
    uint32_t new_index = index ^ ( 1 << target );
    if ( index & ( 1 << target ) )
    {
      new_state.index_to_weight[new_index] += c10[reverse] * weight;
      new_state.index_to_weight[index] += c11[reverse] * weight;
    }
    else
    {
      new_state.index_to_weight[index] += c00[reverse] * weight;
      new_state.index_to_weight[new_index] += c01[reverse] * weight;
    }
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

QRState U2::operator()( const QRState& state, uint32_t target, const bool reverse ) const
{
  QRState new_state;
  for ( const auto& [index, weight] : state.index_to_weight )
  {
    if ( new_state.index_to_weight.find( index ) == new_state.index_to_weight.end() )
      new_state.index_to_weight[index] = 0;
    uint32_t new_index = index ^ ( 1 << target );
    if ( index & ( 1 << target ) )
    {
      new_state.index_to_weight[new_index] += c10[reverse].real() * weight;
      new_state.index_to_weight[index] += c11[reverse].real() * weight;
    }
    else
    {
      new_state.index_to_weight[index] += c00[reverse].real() * weight;
      new_state.index_to_weight[new_index] += c01[reverse].real() * weight;
    }
  }
  for ( auto it = new_state.index_to_weight.begin(); it != new_state.index_to_weight.end(); )
  {
    if ( std::abs( it->second ) < QRState::eps )
      it = new_state.index_to_weight.erase( it );
    else
      ++it;
  }
  new_state.n_bits = state.n_bits;
  return new_state;
}

} // namespace xyz