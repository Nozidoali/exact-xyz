#include "qgate.hpp"
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace xyz
{
std::vector<std::shared_ptr<QGate>> decompose_mcry( const MCRY& gate )
{
  std::vector<std::shared_ptr<QGate>> gates;
  if ( gate.ctrls.size() == 0 )
  {
    gates.push_back( std::make_shared<RY>( gate.target, gate.theta ) );
    return gates;
  }
  if ( gate.ctrls.size() == 1 )
  {
    gates.push_back( std::make_shared<CRY>( gate.ctrls[0], gate.phases[0], gate.theta, gate.target ) );
    return gates;
  }
  return gates;
}

} // namespace xyz