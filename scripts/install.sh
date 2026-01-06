#!/bin/bash
set -e

echo "=== Building XYZ Library and Python Bindings ==="

cd "$(dirname "$0")/.."

if [ ! -d "build" ]; then
    mkdir build
fi

cd build

echo "Configuring CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DBUILD_PYTHON_BINDINGS=ON \
         -DBUILD_XYZ_TESTS=OFF

echo "Building..."
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "✓ Build complete"
echo ""
echo "To use Python bindings:"
echo "  export PYTHONPATH=$(pwd):$(pwd)/../bindings:\$PYTHONPATH"
echo "  python -c 'from xyz_prep import prepare_state'"

