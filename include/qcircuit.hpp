#pragma once
#include "qgate.hpp"
#include "qstate.hpp"
#include <memory>
#include <vector>

namespace xyz
{
class QCircuit
{
private:
  /* data */
  std::vector<std::shared_ptr<QGate>> pGates;
  uint32_t num_qubits = 0;

public:
  // Default constructor
  QCircuit() = default;
  QCircuit( uint32_t num_qubits ) : num_qubits( num_qubits ){};
  void add_gate( std::shared_ptr<QGate> gate );
  int num_cnots() const;
  std::string to_string() const;
};

QCircuit prepare_state( const QState& state );
std::string to_qasm2( const QCircuit& circuit );

} // namespace xyz
