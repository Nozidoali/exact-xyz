#pragma once
#include "qgate.hpp"
#include "qstate.hpp"
#include <memory>
#include <vector>

namespace xyz
{
class QCircuit
{
public:
  uint32_t num_qbits = 0;
  std::vector<std::shared_ptr<QGate>> pGates;

public:
  // Default constructor
  QCircuit() = default;
  QCircuit( uint32_t num_qbits ) : num_qbits( num_qbits ) {};
  void add_gate( std::shared_ptr<QGate> gate );
  uint32_t num_cnots() const;
  std::string to_qasm2() const;
};

QCircuit decompose_circuit( const QCircuit& circuit );

/* State preparation algorithms */
QCircuit prepare_state( const QRState& state );  /* DATE24 */
QCircuit prepare_ghz( uint32_t n );              /* GHZ */

/* ICCAD 24 */
QCircuit resyn( const QCircuit& circuit );

void write_qasm2( const QCircuit& circuit, const std::string& filename );
QCircuit read_qasm2( const std::string& filename, bool verbose = false );


QRState simulate_circuit( const QCircuit& circuit, const QRState& state, bool verbose = false );

} // namespace xyz
