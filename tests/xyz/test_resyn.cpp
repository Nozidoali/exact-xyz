#include "state_test_utils.hpp"

using namespace xyz;
using namespace xyz::testutil;

TEST_CASE("resyn preserves circuit semantics", "[xyz]") {
    QCircuit c(3);
    c.add_gate(std::make_shared<RY>(0, 0.3));
    c.add_gate(std::make_shared<CX>(1, true, 0));
    c.add_gate(std::make_shared<RY>(0, -0.7));
    c.add_gate(std::make_shared<CX>(2, false, 0));
    c.add_gate(std::make_shared<RY>(0, 1.1));

    auto a = simulate_circuit(c, ground_rstate(3), false);
    auto r = resyn(c, false);
    auto b = simulate_circuit(r, ground_rstate(3), false);

    require_close(a, b);
    REQUIRE(r.num_cnots() <= c.num_cnots());
}
