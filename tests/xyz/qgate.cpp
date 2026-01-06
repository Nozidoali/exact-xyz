#include <catch.hpp>
#include <cstdint>
#include "qcircuit.hpp"
#include "qgate.hpp"
#include "qstate.hpp"

using namespace xyz;

TEST_CASE("Test prepare_ghz", "[xyz]") {
    auto circuit = prepare_ghz(4);
    REQUIRE(circuit.num_qbits == 4);
    REQUIRE(circuit.num_cnots() == 3);
    REQUIRE(circuit.lev_cnots() == 3);
}

TEST_CASE("Test prepare_ghz with parallel cnots", "[xyz]") {
    auto circuit = prepare_ghz(4, true);
    REQUIRE(circuit.num_qbits == 4);
    REQUIRE(circuit.num_cnots() == 3);
    REQUIRE(circuit.lev_cnots() == 2);
}

TEST_CASE("Test prepare_w", "[xyz]") {
    auto circuit = prepare_w(4);
    REQUIRE(circuit.num_qbits == 4);
    REQUIRE(circuit.num_cnots() > 0);
    REQUIRE(circuit.lev_cnots() > 0);
}

TEST_CASE("Test prepare_w with logarithmic depth", "[xyz]") {
    auto circuit = prepare_w(4, true);
    REQUIRE(circuit.num_qbits == 4);
    REQUIRE(circuit.num_cnots() > 0);
    REQUIRE(circuit.lev_cnots() > 0);
}