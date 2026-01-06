## exact-xyz

C++17 library for simple quantum circuits, real-amplitude state simulation, and state preparation routines.

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_XYZ_TESTS=ON
cmake --build build -j
```

### Run tests

```bash
./build/tests/run_tests
```

### Tools

```bash
./build/tools/state_simulation -i ../data/input.qasm
./build/tools/prepare_dicke -n 4 -k 2
```

### Python Bindings

Python bindings are available for state preparation functionality. They allow you to generate quantum circuits from numpy arrays and use them with Qiskit.

#### Installation

```bash
# Build with Python bindings enabled
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --build build -j

# Install the Python package
cd bindings
pip install -e .
```

#### Quick Start

```python
import numpy as np
from xyz_prep import prepare_state

# Define a normalized quantum state
coefficients = np.array([0.6, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0])

# Generate the preparation circuit (returns a Qiskit QuantumCircuit)
circuit = prepare_state(coefficients)

print(f"Circuit has {circuit.num_qubits} qubits and {circuit.size()} gates")
print(circuit)
```

See [`bindings/README.md`](bindings/README.md) for detailed documentation and more examples.


