#pragma once

#include <cstdint>
#include <iostream>
#include <map>
#include <optional>
#include <unordered_map>

namespace xyz
{
class QState
{
private:
  /* data */
  std::optional<uint64_t> hash_value = std::nullopt;

public:
  /* data */
  uint32_t n_bits = 0;
  std::map<uint32_t, double> index_to_weight;
  QState() = default;
  QState( std::map<uint32_t, double>& index_to_weight, uint32_t n_bits ) : index_to_weight( index_to_weight ), n_bits( n_bits ){};
  std::unordered_map<uint32_t, double> to_ry_table( uint32_t target ) const;
  bool is_ground() const;
  friend std::ostream& operator<<( std::ostream& os, const QState& obj );
  static constexpr double eps = 1e-6;
  uint64_t repr() const;
  QState clone() const;
};
QState ground_state( uint32_t n_bits );
QState dicke_state( uint32_t n, uint32_t k );
} // namespace xyz
