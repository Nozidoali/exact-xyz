#include "synthesis.hpp"

namespace xyz
{
Schedule layout_synthesis( const QCircuit& circuit, const Config& config )
{
  // we first transpile the circuit to a {CNOT, U2} basis
  QCircuit transpiled_circuit = decompose_circuit( circuit );

  Schedule schedule;
  return schedule;
}
} // namespace xyz