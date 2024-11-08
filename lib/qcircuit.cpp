#include "qcircuit.hpp"
#include <fstream>
#include <memory>
#include <string>
#include <typeinfo>

namespace xyz
{
void QCircuit::add_gate( std::shared_ptr<QGate> gate )
{
  pGates.push_back( gate );
}
std::string QCircuit::to_qasm2() const
{
  std::string qasm = "";
  qasm += "OPENQASM 2.0;\n";
  qasm += "include \"qelib1.inc\";\n";
  qasm += "qreg q[" + std::to_string( num_qbits ) + "];\n";
  for ( const auto& pGate : pGates )
  {
    qasm += pGate->to_string() + ";\n";
  }
  return qasm;
}
QCircuit decompose_circuit( const QCircuit& circuit )
{
  QCircuit new_circuit( circuit.num_qbits );
  for ( const auto& pGate : circuit.pGates )
  {
    // if the class is a cry
    std::shared_ptr<CRY> cry_gate = std::dynamic_pointer_cast<CRY>( pGate );
    if ( cry_gate )
    {
      new_circuit.add_gate( std::make_shared<RY>( cry_gate->target, cry_gate->theta / 2 ) );
      new_circuit.add_gate( std::make_shared<CX>( cry_gate->ctrl, true, cry_gate->target ) );
      if ( cry_gate->phase == false )
        new_circuit.add_gate( std::make_shared<RY>( cry_gate->target, cry_gate->theta / 2 ) );
      else
        new_circuit.add_gate( std::make_shared<RY>( cry_gate->target, -cry_gate->theta / 2 ) );
      new_circuit.add_gate( std::make_shared<CX>( cry_gate->ctrl, true, cry_gate->target ) );
      continue;
    }
    // if the class is a cx
    std::shared_ptr<CX> cx_gate = std::dynamic_pointer_cast<CX>( pGate );
    if ( cx_gate )
    {
      if ( cx_gate->phase == false )
      {
        new_circuit.add_gate( std::make_shared<X>( cx_gate->ctrl ) );
        new_circuit.add_gate( std::make_shared<CX>( cx_gate->ctrl, true, cx_gate->target ) );
        new_circuit.add_gate( std::make_shared<X>( cx_gate->ctrl ) );
        continue;
      }
      else
      {
        new_circuit.add_gate( pGate );
        continue;
      }
    }

    new_circuit.add_gate( pGate );
  }
  return new_circuit;
}
uint32_t QCircuit::num_cnots() const
{
  int count = 0;
  for ( const auto& pGate : pGates )
  {
    count += pGate->num_cnots();
  }
  return count;
}
void write_qasm2( const QCircuit& circuit, const std::string& filename )
{
  std::ofstream file;
  file.open( filename );
  file << circuit.to_qasm2();
  file.close();
}
QCircuit read_qasm2( const std::string& filename, bool verbose )
{
  std::ifstream file( filename );
  std::string line;
  QCircuit circuit;

  if ( verbose )
    std::cout << "Reading file: " << filename << std::endl;

  while ( std::getline( file, line ) )
  {
    if ( verbose )
    {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << line << std::endl;
    }
    if ( line.find( "qreg " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto num_qbits = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      circuit.num_qbits = num_qbits;
      if ( verbose )
        std::cout << "num_qbits: " << num_qbits << std::endl;
      continue;
    }
    else if ( line.find( "x " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto target = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      circuit.add_gate( std::make_shared<X>( target ) );
      continue;
    }
    else if ( line.find( "cx " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "[", pos2 );
      auto pos4 = line.find( "]", pos3 );
      auto ctrl = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto target = std::stoi( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      circuit.add_gate( std::make_shared<CX>( ctrl, true, target ) );
      continue;
    }
    else if ( line.find( "ry" ) == 0 ) // ry
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "(" );
      auto pos4 = line.find( ")", pos3 );
      auto target = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto theta = std::stod( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      circuit.add_gate( std::make_shared<RY>( target, theta ) );
      continue;
    }
    else if ( line.find( "cry" ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "[", pos2 );
      auto pos4 = line.find( "]", pos3 );
      auto pos5 = line.find( "(" );
      auto pos6 = line.find( ")", pos5 );
      auto ctrl = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto target = std::stoi( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      auto theta = std::stod( line.substr( pos5 + 1, pos6 - pos5 - 1 ) );
      circuit.add_gate( std::make_shared<CRY>( ctrl, true, theta, target ) );
      continue;
    }
    else if ( line.find( "h " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto target = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      circuit.add_gate( std::make_shared<H>( target ) );
      continue;
    }
    else if ( line.find( "t " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto target = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      circuit.add_gate( std::make_shared<T>( target ) );
      continue;
    }
    else if ( line.find( "cx_false" ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "[", pos2 );
      auto pos4 = line.find( "]", pos3 );
      auto ctrl = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto target = std::stoi( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      circuit.add_gate( std::make_shared<CX>( ctrl, false, target ) );
      continue;
    }
    else if ( line.find( "cry_false" ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "[", pos2 );
      auto pos4 = line.find( "]", pos3 );
      auto pos5 = line.find( "(" );
      auto pos6 = line.find( ")", pos5 );
      auto ctrl = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto target = std::stoi( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      auto theta = std::stod( line.substr( pos5 + 1, pos6 - pos5 - 1 ) );
      circuit.add_gate( std::make_shared<CRY>( ctrl, false, theta, target ) );
      continue;
    }
    else if ( line.find( "ccx " ) == 0 )
    {
      auto pos1 = line.find( "[" );
      auto pos2 = line.find( "]", pos1 );
      auto pos3 = line.find( "[", pos2 );
      auto pos4 = line.find( "]", pos3 );
      auto pos5 = line.find( "[", pos4 );
      auto pos6 = line.find( "]", pos5 );
      auto ctrl1 = std::stoi( line.substr( pos1 + 1, pos2 - pos1 - 1 ) );
      auto ctrl2 = std::stoi( line.substr( pos3 + 1, pos4 - pos3 - 1 ) );
      auto target = std::stoi( line.substr( pos5 + 1, pos6 - pos5 - 1 ) );
      circuit.add_gate( std::make_shared<CCX>( ctrl1, ctrl2, target ) );
      continue;
    }
    else
    {
      std::cerr << "Unknown gate: " << line << std::endl;
    }
  }
  return circuit;
}
QRState simulate_circuit( const QCircuit& circuit, const QRState& state, bool verbose )
{
  QRState new_state = state;
  for ( const auto& pGate : circuit.pGates )
  {
    if ( verbose )
    {
      std::cout << "Applying gate: " << pGate->to_string() << std::endl;
      std::cout << "State before: " << new_state << std::endl;
    }
    new_state = ( *pGate )( new_state );
    if ( verbose )
    {
      std::cout << "State after: " << new_state << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
  }
  return new_state;
}
} // namespace xyz
