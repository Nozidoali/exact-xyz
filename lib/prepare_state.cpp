#include "qcircuit.hpp"
#include <algorithm>
#include <cassert>
#include <fmt/format.h>
#include <iostream>
#include <memory>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace xyz
{
namespace detail
{
struct bfs_state
{
  uint32_t mem_idx;
  uint32_t cnot_cost;
  bool operator<( const bfs_state& other ) const
  {
    // TODO: A* search
    return cnot_cost > other.cnot_cost;
  }
  bfs_state( uint32_t mem_idx, uint32_t cnot_cost ) : mem_idx( mem_idx ), cnot_cost( cnot_cost ){};
};

struct memorized_state
{
  uint32_t prev;
  QRState state;
  std::shared_ptr<QGate> pGate = nullptr;
  memorized_state( const QRState& state, uint32_t prev ) : state( state ), prev( prev ){};
  memorized_state( const QRState& state, uint32_t prev, std::shared_ptr<QGate> pGate ) : state( state ), prev( prev ), pGate( pGate ){};
};

std::vector<std::shared_ptr<QGate>> enumerate_gates( const QRState& state )
{
  std::vector<std::shared_ptr<QGate>> gates;
  if ( state.index_to_weight.size() == 1 )
  {
    auto index = state.index_to_weight.begin()->first;
    for ( uint32_t target = 0; target < state.n_bits; target++ )
    {
      if ( ( index >> target ) & 1 )
      {
        gates.push_back( std::make_shared<X>( target ) );
        return gates;
      }
    }
  }
  for ( uint32_t target = 0; target < state.n_bits; target++ )
  {
    const auto ry_table = state.to_ry_table( target );
    std::optional<double> theta;
    for ( const auto& [index, _theta] : ry_table )
    {
      if ( theta.has_value() && theta.value() != _theta )
      {
        theta.reset();
        break;
      }
      theta = _theta;
    }
    if ( theta.has_value() && !Rotation::is_trivial( theta.value(), true ) )
    {
      gates.push_back( std::make_shared<RY>( target, theta.value() ) );
      return gates;
    }
  }
  for ( uint32_t target = 0; target < state.n_bits; target++ )
  {
    const auto ry_table = state.to_ry_table( target );
    for ( uint32_t ctrl = 0; ctrl < state.n_bits; ctrl++ )
      for ( bool phase : { false, true } )
      {
        if ( ctrl == target )
          continue;
        std::optional<double> theta = std::nullopt;
        for ( const auto& [index, _theta] : ry_table )
        {
          if ( (bool)( ( index >> ctrl ) & 1 ) == phase )
          {
            if ( theta.has_value() && theta.value() != _theta )
            {
              theta.reset();
              break;
            }
            theta = _theta;
          }
        }
        if ( theta.has_value() && !Rotation::is_trivial( theta.value(), true ) )
        {
          gates.push_back( std::make_shared<CRY>( ctrl, phase, theta.value(), target ) );
          gates.push_back( std::make_shared<CRY>( ctrl, phase, -M_PI + theta.value(), target ) );
        }
      }
  }
  for ( uint32_t target = 0; target < state.n_bits; target++ )
    for ( uint32_t ctrl = 0; ctrl < state.n_bits; ctrl++ )
    {
      if ( ctrl == target )
        continue;
      gates.push_back( std::make_shared<CX>( ctrl, true, target ) );
    }
  return gates;
}

template<bool verbose>
void prepare_state_impl( const QRState& state, QCircuit& circuit )
{
  std::priority_queue<bfs_state> q;
  std::vector<memorized_state> states;
  std::unordered_map<uint64_t, uint32_t> state_to_cost;
  std::unordered_set<uint64_t> visited;
  states.push_back( memorized_state( state, -1 ) );
  q.push( bfs_state( 0, 0 ) );
  state_to_cost[state.repr()] = 0;
  std::optional<uint32_t> solution = std::nullopt;
  while ( !q.empty() )
  {
    auto e = q.top();
    q.pop();
    auto current_state = states[e.mem_idx];
    if ( visited.find( current_state.state.repr() ) != visited.end() )
      continue;
    visited.insert( current_state.state.repr() );
    if ( current_state.state.is_ground() )
    {
      solution = e.mem_idx;
      break;
    }
    if constexpr ( verbose )
      fmt::print( "current_state: {} repr = {}\n", current_state.state.to_string(), current_state.state.repr() );
    for ( const auto& gate : enumerate_gates( current_state.state ) )
    {
      if constexpr ( verbose )
        std::cout << "gate: " << *gate << std::endl;
      auto new_state = ( *gate )( current_state.state, true );
      auto new_cost = e.cnot_cost + gate->num_cnots();
      auto new_repr = new_state.repr();
      if ( visited.find( new_repr ) != visited.end() )
        continue;
      if ( state_to_cost.find( new_repr ) == state_to_cost.end() || state_to_cost[new_repr] > e.cnot_cost + 1 )
      {
        if constexpr ( verbose )
          std::cout << "\tadd: " << *gate << " -> " << new_state << std::endl;
        state_to_cost[new_repr] = new_cost;
        q.push( bfs_state( states.size(), new_cost ) );
        states.push_back( memorized_state( new_state, e.mem_idx, gate ) );
      }
    }
  }
  if ( !solution.has_value() )
  {
    std::cerr << "No solution found!" << std::endl;
    return;
  }
  // back trace
  for ( uint32_t idx = solution.value(); idx != 0; idx = states[idx].prev )
  {
    if constexpr ( verbose )
      std::cout << "gate: " << *states[idx].pGate << " -> " << states[idx].state << std::endl;
    circuit.add_gate( states[idx].pGate );
  }
}

template<bool verbose>
std::pair<uint32_t, uint32_t> maximize_difference_once( uint32_t num_qubits, std::unordered_set<uint32_t>& indices, std::unordered_map<uint32_t, bool>& diff_values )
{
  int max_diff = -1;
  std::unordered_set<uint32_t> max_diff_indices_1;
  uint32_t max_diff_qubit = 0;
  uint32_t max_diff_value = false;
  uint32_t length = indices.size();
  std::unordered_set<uint32_t> indices_1;

  for ( uint32_t qubit = 0; qubit < num_qubits; qubit++ )
  {
    if ( diff_values.find( qubit ) != diff_values.end() )
      continue;
    indices_1.clear();
    for ( auto index : indices )
      if ( ( index >> qubit ) & (uint32_t)1 )
        indices_1.insert( index );

    int diff = abs( (int)length - 2 * (int)indices_1.size() );
    if ( diff == length )
      continue;
    if constexpr ( verbose )
      fmt::print( "qubit: {}, diff: {}\n", qubit, diff );

    if ( diff > max_diff )
    {
      max_diff = diff;
      max_diff_indices_1 = indices_1;
      max_diff_qubit = qubit;
      max_diff_value = length > 2 * indices_1.size();
    }
    if ( max_diff == length - 1 )
      break;
  }

  if ( max_diff_value )
    indices = max_diff_indices_1;
  else
  {
    // indices = indices - indices_1
    for ( auto index : max_diff_indices_1 )
      indices.erase( index );
  }
  diff_values[max_diff_qubit] = max_diff_value;
  return { max_diff_qubit, max_diff_value };
}

template<bool verbose>
QRState cardinality_reduction_impl( QCircuit& qcircuit, const QRState& state )
{
  if constexpr ( verbose )
  {
    fmt::print( "----------------------------------------\n" );
    fmt::print( "[i] state: {}\n", state.to_string() );
    fmt::print( "\tcardinality: {}\n", state.cardinality() );
  }
  QRState new_state = state.clone();
  std::unordered_set<uint32_t> indices;
  for ( auto [index, weight] : state.index_to_weight )
    indices.insert( index );
  std::unordered_map<uint32_t, bool> diff_values;
  uint32_t diff_qubit, diff_value;
  while ( indices.size() > 1 )
  {
    std::tie( diff_qubit, diff_value ) = maximize_difference_once<verbose>( state.n_bits, indices, diff_values );
    if constexpr ( verbose )
    {
      fmt::print( "[done] diff_qubit: {}, diff_value: {}, indices.size(): {}\n", diff_qubit, diff_value, indices.size() );
    }
  }
  uint32_t index0 = *indices.begin();
  diff_values.erase( diff_qubit );

  if constexpr ( verbose )
  {
    fmt::print( "[i] found index0 = {}\n", index0 );
    fmt::print( "diff_values: " );
    for ( auto [qubit, value] : diff_values )
      fmt::print( "{}: {}, ", qubit, value );
    fmt::print( "\n" );
  }
  std::unordered_set<uint32_t> candidates;
  for ( auto [index, weight] : state.index_to_weight )
  {
    if ( indices.find( index ) != indices.end() )
      continue;
    bool valid = true;
    for ( auto [qubit, value] : diff_values )
      if ( ( ( index >> qubit ) & (uint32_t)1 ) != value )
      {
        valid = false;
        break;
      }
    if constexpr ( verbose )
      fmt::print( "index: {}, valid: {}\n", index, valid );
    if ( valid )
      candidates.insert( index );
  }
  if constexpr ( verbose )
  {
    fmt::print( "[i] candidates for index1: " );
    for ( auto index : candidates )
      fmt::print( "{}, ", index );
    fmt::print( "\n" );
  }
  while ( candidates.size() > 1 )
    (void)maximize_difference_once<verbose>( state.n_bits, candidates, diff_values );
  uint32_t index1 = *candidates.begin();
  for ( uint32_t qubit = 0; qubit < state.n_bits; qubit++ )
  {
    if ( ( index0 >> qubit & 1 ) == ( index1 >> qubit & 1 ) )
      continue;
    if ( qubit == diff_qubit )
      continue;
    qcircuit.add_gate( std::make_shared<CX>( diff_qubit, diff_value, qubit ) );
    new_state = ( *std::make_shared<CX>( diff_qubit, diff_value, qubit ) )( new_state, true );
  }

  std::vector<uint32_t> ctrls;
  std::vector<bool> phases;
  for ( auto [qubit, value] : diff_values )
  {
    ctrls.push_back( qubit );
    phases.push_back( value );
  }
  index0 = index1 ^ ( (uint32_t)1 << diff_qubit );
  if constexpr ( verbose )
  {
    fmt::print( "[i] new state before merging: {}\n", new_state.to_string() );
    fmt::print( "\tdiff_qubit: {}, diff_value: {}\n", diff_qubit, diff_value );
    fmt::print( "[i] merging indices from index0 = {} to index1 = {}\n", index0, index1 );
  }
  uint32_t idx0 = index1 & ( ~( (uint32_t)1 << diff_qubit ) );
  uint32_t idx1 = index1 | ( ( (uint32_t)1 << diff_qubit ) );
  auto it0 = new_state.index_to_weight.find( idx0 );
  auto it1 = new_state.index_to_weight.find( idx1 );
  assert( it0 != new_state.index_to_weight.end() );
  assert( it1 != new_state.index_to_weight.end() );
  double theta = 2.0 * atan2l( (long double)it1->second, (long double)it0->second );
  if constexpr ( verbose )
  {
    fmt::print( "[i] weight0: {}, weight1: {}, theta: {}\n", it0->second, it1->second, theta );
    for ( auto [qubit, value] : diff_values )
      fmt::print( " qubit = {}, phase = {} \n", qubit, value );
  }
  if ( index1 & ( 1 << diff_qubit ) )
    theta = -M_PI + theta;
  assert( theta < 2 * M_PI );
  qcircuit.add_gate( std::make_shared<MCRY>( ctrls, phases, theta, diff_qubit ) );
  new_state = ( *std::make_shared<MCRY>( ctrls, phases, theta, diff_qubit ) )( new_state, true );
  if constexpr ( verbose )
    fmt::print( "[i] new state after merging: {}\n", new_state.to_string() );
  return new_state;
}

// Adds a CX gate to the circuit, modifying qVal and counting the CNOTs
void add_cx( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control, uint32_t target )
{
  if ( qVal[control].has_value() && !qVal[control].value() )
    return;
  if ( qVal[control].has_value() && qVal[control].value() )
  {
    circuit.add_gate( std::make_shared<X>( target ) );
    if ( qVal[target].has_value() )
      qVal[target] = !qVal[target].value();
  }
  else
  {
    circuit.add_gate( std::make_shared<CX>( control, true, target ) );
    qVal[target] = std::nullopt;
  }
}

// Adds an RY rotation gate to the circuit
void add_ry( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t target, double theta )
{
  circuit.add_gate( std::make_shared<RY>( target, theta ) );
  qVal[target] = std::nullopt;
}

// Adds a conditional RY (CRY) rotation gate with a control qubit
void add_cry( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control, uint32_t target, double theta )
{
  if ( qVal[control].has_value() && !qVal[control].value() )
    return;
  if ( qVal[control].has_value() && qVal[control].value() )
  {
    add_ry( circuit, qVal, target, theta );
  }
  else
  {
    circuit.add_gate( std::make_shared<CRY>( control, true, theta, target ) );
    qVal[target] = std::nullopt;
  }
}

// Adds a multi-controlled RY (MCRY) rotation gate with two control qubits
void add_mcry( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t control1, uint32_t control2, uint32_t target, double theta )
{
  if ( qVal[control1].has_value() && !qVal[control1].value() )
    return;
  if ( qVal[control2].has_value() && !qVal[control2].value() )
    return;
  if ( qVal[control1].has_value() && qVal[control1].value() )
  {
    add_cry( circuit, qVal, control2, target, theta );
    return;
  }
  if ( qVal[control2].has_value() && qVal[control2].value() )
  {
    add_cry( circuit, qVal, control1, target, theta );
    return;
  }
  std::vector<uint32_t> controls = { control1, control2 };
  circuit.add_gate( std::make_shared<MCRY>( controls, theta, target ) );
  qVal[target] = std::nullopt;
}

// Function to insert a "μ" rotation gate in the circuit
void insert_mu( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t k, uint32_t j )
{
  double theta = 2 * std::acos( std::sqrt( 1.0 / ( n - j ) ) );
  add_cx( circuit, qVal, j + 1, j );
  add_cry( circuit, qVal, j, j + 1, theta );
  add_cx( circuit, qVal, j + 1, j );
}

// Function to insert an "M" rotation gate in the circuit
void insert_M( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t k, uint32_t j, uint32_t i )
{
  double theta = 2 * std::acos( std::sqrt( ( i + 1.0 ) / ( n - j ) ) );
  add_cx( circuit, qVal, j + i + 1, j );
  add_mcry( circuit, qVal, j + i, j, j + i + 1, theta );
  add_cx( circuit, qVal, j + i + 1, j );
}

// Function to insert a single-control controlled sequence (scs) into the circuit
void insert_scs( QCircuit& circuit, std::vector<std::optional<bool>>& qVal, uint32_t n, uint32_t k, uint32_t j )
{
  assert( k >= 1 && "k should be greater than or equal to 1" );
  insert_mu( circuit, qVal, n, k, j );
  for ( uint32_t i = 1; i < k; ++i )
  {
    if ( j + i + 1 >= n )
      break;
    insert_M( circuit, qVal, n, k, j, i );
  }
}

} // namespace detail

QCircuit prepare_state( const QRState& state )
{
  QCircuit circuit( state.n_bits );
  detail::prepare_state_impl<true>( state, circuit );
  return circuit;
}

QCircuit prepare_ghz( uint32_t n, bool log_depth )
{
  /* Reference: https://arxiv.org/abs/1807.05572 */
  QCircuit circuit( n );
  circuit.add_gate( std::make_shared<H>( 0 ) );
  if ( log_depth ) /* logarithmic depth */
    for ( uint32_t i = 1u; i < n; i <<= 1 )
      for ( uint32_t j = 0; j < i && j + i < n; j++ )
        circuit.add_gate( std::make_shared<CX>( j, true, j + i ) );
  else /* linear depth */
    for ( uint32_t i = 1; i < n; i++ )
      circuit.add_gate( std::make_shared<CX>( 0, true, i ) );
  return circuit;
}

QCircuit prepare_w( uint32_t n, bool log_depth, bool cnot_opt )
{
  /* Reference: https://arxiv.org/abs/1807.05572 */
  QCircuit circuit( n );
  circuit.add_gate( std::make_shared<X>( 0 ) );
  if ( !log_depth ) /* linear depth */
  {
    for ( uint32_t i = 1; i < n; i++ )
    {
      uint32_t j = i - 1;
      double p = 1.0 / (double)( n - j );
      double theta = 2 * std::atan2( std::sqrt( 1 - p ), std::sqrt( p ) );
      if ( cnot_opt )
      {
        // replace it using two RYs and a CNOT
        circuit.add_gate( std::make_shared<RY>( j, -( theta - M_PI ) / 2 ) );
        circuit.add_gate( std::make_shared<CX>( j, true, i ) );
        circuit.add_gate( std::make_shared<RY>( j, ( theta - M_PI ) / 2 ) );
      }
      else
        circuit.add_gate( std::make_shared<CRY>( j, true, theta, i ) );
      circuit.add_gate( std::make_shared<CX>( i, true, j ) );
    }
  }
  else /* logarithmic depth */
  {
    using P = std::pair<uint32_t, uint32_t>;
    std::queue<std::array<uint32_t, 3>> dicotomies;
    dicotomies.push( { { 0, n, n >> 1 } } ); /* q, tot, curr */
    uint32_t q_next = 1;
    while ( !dicotomies.empty() )
    {
      auto [q, total, curr] = dicotomies.front();
      dicotomies.pop();
      if ( total < 2 )
        continue;
      uint32_t total_l = total >> 1;
      uint32_t curr_l = curr >> 1; /* l <= r */
      uint32_t total_r = total - total_l;
      uint32_t curr_r = curr - curr_l;
      if ( total_l == 1 && curr_l == 1 ) /* the only corner case */
        dicotomies.push( { q, total_r, curr_r } );
      else /* divide and conquer */
      {
        dicotomies.push( { q, total_l, curr_l } );
        dicotomies.push( { q_next, total_r, curr_r } );
      }
      double p = (double)curr / (double)total;
      double theta = 2 * std::atan2( std::sqrt( 1 - p ), std::sqrt( p ) );
      if ( cnot_opt )
      {
        // replace it using two RYs and a CNOT
        circuit.add_gate( std::make_shared<RY>( q, -( theta - M_PI ) / 2 ) );
        circuit.add_gate( std::make_shared<CX>( q, true, q_next ) );
        circuit.add_gate( std::make_shared<RY>( q, ( theta - M_PI ) / 2 ) );
      }
      else
        circuit.add_gate( std::make_shared<CRY>( q, true, theta, q_next ) );
      circuit.add_gate( std::make_shared<CX>( q_next, true, q ) );
      q_next++;
    }
  }
  return circuit;
}

QCircuit prepare_sparse_state( const QRState& state )
{
  QCircuit circuit( state.n_bits );
  QRState curr_state = state.clone();
  while ( curr_state.cardinality() > 1 )
    curr_state = detail::cardinality_reduction_impl<true>( circuit, curr_state );
  uint32_t index = curr_state.index_to_weight.begin()->first;
  for ( uint32_t qubit = 0; qubit < state.n_bits; qubit++ )
    if ( ( index >> qubit ) & (uint32_t)1 )
      circuit.add_gate( std::make_shared<X>( qubit ) );
  circuit.reverse();
  return circuit;
}

QCircuit prepare_dicke_state( int n, int k )
{
  QCircuit circuit( n );
  std::vector<std::optional<bool>> qVal( n, false );
  for ( int i = 0; i < k; ++i )
  {
    circuit.add_gate( std::make_shared<X>( i ) );
    qVal[i] = true;
  }
  for ( int i = 0; i < n - 1; ++i )
    detail::insert_scs( circuit, qVal, n, k, i );
  return circuit;
}

} // namespace xyz