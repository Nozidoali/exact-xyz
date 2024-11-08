#include <qcircuit.hpp>

using namespace xyz;

int main()
{
  auto qc = read_qasm2( "../data/qasm/tof_3.qasm", true );
  std::cout << "Number of CNOTs: " << qc.num_cnots() << std::endl;
  write_qasm2( qc, "output.qasm" );
  return 0;
}