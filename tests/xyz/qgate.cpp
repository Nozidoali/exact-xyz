#include <catch.hpp>
#include <xyz.hpp>

using namespace xyz;

TEST_CASE( "Test prepare_ghz", "[xyz]" )
{
  auto circuit = prepare_ghz( 3 );
  REQUIRE( circuit.num_qbits == 3 );
  REQUIRE( circuit.num_cnots() == 2 );
}