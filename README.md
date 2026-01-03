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


