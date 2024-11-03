#include "qcircuit.hpp"
#include <algorithm>
#include <optional>
namespace xyz
{
namespace detail
{
using RLUT = std::vector<std::pair<uint32_t, std::pair<double, double>>>;
using CXT = std::vector<std::pair<int, bool>>;
using RSOL = std::vector<double>;
void enumerate_cnot_templates( std::vector<CXT>& templates, const std::vector<uint32_t>& controls, uint32_t n_cnots, uint32_t n_cnots_added, CXT& curr_template )
{
  if ( n_cnots_added == n_cnots )
  {
    templates.push_back( curr_template );
    return;
  }
  for ( auto control : controls )
    for ( auto phase : { true, false } )
    {
      curr_template.push_back( std::make_pair( control, phase ) );
      enumerate_cnot_templates( templates, controls, n_cnots, n_cnots_added + 1, curr_template );
      curr_template.pop_back();
    }
}
std::optional<RSOL> gaussian_elimination( std::vector<std::vector<int>>& R, std::vector<double>& b )
{
  uint32_t m = R.size();    /* number of rows */
  uint32_t n = R[0].size(); /* number of columns */

  // print the R matrix and the RHS
  for ( uint32_t i = 0; i < m; i++ )
  {
    for ( uint32_t j = 0; j < n; j++ )
      std::cout << R[i][j] << " ";
    std::cout << " | " << b[i] << std::endl;
  }
  std::cout << std::endl;

  for ( uint32_t i = 0; i < m; i++ )
  {
    /* find the pivot */
    uint32_t pivot = i;
    for ( uint32_t j = i; j < m; j++ )
      if ( R[j][i] != 0 )
      {
        pivot = j;
        break;
      }
    if ( pivot != i )
    {
      std::swap( R[i], R[pivot] );
      std::swap( b[i], b[pivot] );
    }
    if ( R[i][i] == 0 )
      return std::nullopt;
    /* eliminate the column */
    for ( uint32_t j = i + 1; j < m; j++ )
    {
      if ( R[j][i] == 0 )
        continue;
      double ratio = (double)R[j][i] / R[i][i];
      for ( uint32_t k = i; k < n; k++ )
        R[j][k] -= R[i][k] * ratio;
      b[j] -= b[i] * ratio;
    }
  }

  RSOL sol;
  for ( int i = m - 1; i >= 0; i-- )
  {
    double sum = 0;
    for ( uint32_t j = i + 1; j < n; j++ )
      sum += R[i][j] * sol[j - i - 1];
    sol.push_back( ( b[i] - sum ) / R[i][i] );
  }

  std::cout << "solution: ";
  for ( auto elem : sol )
    std::cout << elem << " ";
  std::cout << std::endl;

  return sol;
}

template<bool verbose>
std::optional<std::pair<RSOL, CXT>> rotation_solver( const RLUT& rlut, const std::vector<uint32_t>& controls, uint32_t max_cnots )
{
  RSOL sol;
  std::vector<CXT> templates;
  for ( uint32_t n_cnots = 0; n_cnots < max_cnots; ++n_cnots )
  {
    templates.clear();
    CXT _config;
    enumerate_cnot_templates( templates, controls, n_cnots, 0, _config );
    // calculate the Polar Matrix of each template, given the indices in RLUT
    uint32_t m = rlut.size();       /* cardinality of the care set */
    uint32_t k = n_cnots + 1;       /* number of rotations in RSOL */
    for ( auto config : templates ) /* enumerate template */
    {
      if constexpr ( verbose )
      {
        for ( uint32_t i = 0; i < 80; i++ )
          std::cout << "-";
        std::cout << std::endl;
        std::cout << "config: ";
        for ( auto elem : config )
          std::cout << "q" << elem.first << " " << elem.second << " ";
        std::cout << std::endl;
      }
      std::vector<std::vector<int>> R; /* Polar Matrix */
      std::vector<double> b;           /* RHS */
      for ( uint32_t i = 0; i < m; i++ )
      {
        uint32_t index = rlut[i].first;
        if constexpr ( verbose )
          std::cout << "index: " << index << std::endl;
        int polarity = 1;
        std::vector<int> row;
        row.push_back( polarity );
        for ( int j = k - 2; j >= 0; --j )
        {
          uint32_t control = config[j].first;
          uint32_t phase = (uint32_t)config[j].second;
          if constexpr ( verbose )
            std::cout << "j: " << j << " control: " << control << " phase: " << phase << " _phase: " << ( ( index >> control ) & 1u ) << std::endl;
          if ( ( ( index >> control ) & 1u ) == phase )
            polarity *= -1;
          row.push_back( polarity );
        }
        double coeff;
        coeff = ( polarity == 1 ) ? 0 : M_PI;
        coeff += rlut[i].second.first * polarity;
        coeff -= rlut[i].second.second;
        b.push_back( coeff );
        std::reverse( row.begin(), row.end() );
        R.push_back( row );
      }
      if constexpr ( verbose ) /* print the R matrix */
      {
        for ( auto row : R )
        {
          for ( auto elem : row )
            std::cout << elem << " ";
          std::cout << std::endl;
        }
      }
      auto _sol = gaussian_elimination( R, b );
      if ( _sol )
        return std::make_pair( *_sol, config );
    }
  }
  return std::nullopt;
}

template<bool verbose>
void resyn_impl( const QCircuit& circuit, QCircuit& new_circuit, const QState& initial_state, uint32_t pos )
{
  if ( pos >= circuit.pGates.size() )
    return;
  if constexpr ( verbose )
    std::cout << "initial_state: " << initial_state << std::endl;
  auto gate = circuit.pGates[pos];
  uint32_t target = gate->target;
  uint32_t initial_cost = gate->num_cnots();
  QState new_state = ( *gate )( initial_state );
  if constexpr ( verbose )
    std::cout << "gate: " << *gate << std::endl;
  uint32_t new_pos = pos + 1;
  for ( ; new_pos < circuit.pGates.size(); new_pos++ )
  {
    auto _gate = circuit.pGates[new_pos];
    if ( _gate->target != target )
      break;
    initial_cost += _gate->num_cnots();
    new_state = ( *_gate )( new_state );
    if constexpr ( verbose )
      std::cout << "gate: " << *_gate << std::endl;
  }
  if constexpr ( verbose )
    std::cout << "final_state: " << new_state << std::endl;

  auto initial_ry = initial_state.to_ry_table( target );
  auto final_ry = new_state.to_ry_table( target );
  RLUT rotation_lut = RLUT();
  for ( auto it : initial_ry )
  {
    auto _it = final_ry.find( it.first );
    rotation_lut.push_back( std::make_pair( it.first, std::make_pair( it.second, _it->second ) ) );
    if constexpr ( verbose )
      std::cout << it.first << " " << it.second << " " << _it->second << std::endl;
  }
  if constexpr ( verbose )
    std::cout << std::endl;

  /* run resynthesis and try resubstitution */
  std::vector<uint32_t> controls;
  for ( uint32_t i = 0; i < circuit.num_qbits; i++ )
    if ( i != target )
      controls.push_back( i );
  auto sol = rotation_solver<true>( rotation_lut, controls, initial_cost );
  if ( sol )
  {
    std::cout << "solution found" << std::endl;

    RSOL rotation_angles = sol->first;
    CXT config = sol->second;
    uint32_t n_cnots = config.size();
    std::cout << "n_cnots: " << n_cnots << std::endl;
    for ( uint32_t i = 0; i < n_cnots; i++ )
    {
      auto control = config[i].first;
      auto phase = config[i].second;
      new_circuit.add_gate( std::make_shared<RY>( target, rotation_angles[i] ) );
      new_circuit.add_gate( std::make_shared<CX>( control, phase, target ) );
    }
    new_circuit.add_gate( std::make_shared<RY>( target, rotation_angles[n_cnots] ) );
  }
  else
  {
    for ( uint32_t i = pos; i < new_pos; i++ )
      new_circuit.add_gate( circuit.pGates[i] );
  }
  resyn_impl<verbose>( circuit, new_circuit, new_state, new_pos );
}
} // namespace detail

QCircuit resyn( const QCircuit& circuit )
{
  QCircuit new_circuit( circuit.num_qbits );
  detail::resyn_impl<false>( circuit, new_circuit, ground_state( circuit.num_qbits ), 0 );
  return new_circuit;
}

} // namespace xyz