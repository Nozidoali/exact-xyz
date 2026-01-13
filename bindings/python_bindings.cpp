#include "prepare-state.hpp"
#include "qcircuit.hpp"
#include "qstate.hpp"
#include "transpile.hpp"

#include <cmath>
#include <map>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

std::string prepare_state_from_array(py::array_t<double> coefficients, double eps = 1e-3, bool verbose = false) {
    py::buffer_info buf = coefficients.request();

    if (buf.ndim != 1)
        throw std::invalid_argument("Input must be 1D array");

    size_t length = buf.shape[0];
    if (length == 0 || (length & (length - 1)) != 0)
        throw std::invalid_argument("Length must be power of 2");

    uint32_t n_bits = 0;
    size_t   temp   = length;
    while (temp > 1) {
        temp >>= 1;
        n_bits++;
    }

    if (n_bits >= 32)
        throw std::invalid_argument("Too many qubits (max 31)");

    double* ptr          = static_cast<double*>(buf.ptr);
    double  norm_squared = 0.0;
    for (size_t i = 0; i < length; i++)
        norm_squared += ptr[i] * ptr[i];

    if (std::abs(norm_squared - 1.0) > 1e-4)
        throw std::invalid_argument("State not normalized");

    std::map<uint32_t, double> index_to_weight;
    for (size_t i = 0; i < length; i++)
        if (std::abs(ptr[i]) >= xyz::QRState::eps)
            index_to_weight[static_cast<uint32_t>(i)] = ptr[i];

    if (index_to_weight.empty())
        throw std::invalid_argument("All coefficients too small");

    xyz::QRState  state(index_to_weight, n_bits);
    xyz::QCircuit circuit    = xyz::prepare_state_auto(state, verbose);
    xyz::QCircuit decomposed = xyz::decompose_circuit(circuit);
    xyz::QCircuit transpiled = xyz::transpile_clifford_t(decomposed, eps);

    return transpiled.to_qasm2();
}

PYBIND11_MODULE(xyz_bindings, m) {
    m.def("prepare_state_from_array", &prepare_state_from_array, py::arg("coefficients"), py::arg("eps") = 1e-3,
          py::arg("verbose") = false);
    m.attr("__version__") = "0.1.0";
}
