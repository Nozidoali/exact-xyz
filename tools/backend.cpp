#include <arch.hpp>
#include <cstdint>
#include <iostream>
#include <map>
#include <qcircuit.hpp>
#include <arch.hpp>

using namespace xyz;
using namespace xyz_backend;
int main() {
    Config config = parseConfig("../data/arch.json");

    // Output parsed data for verification
    std::cout << "Configuration Name: " << config.name << "\n";
    std::cout << "Number of Storage Zones: " << config.storage_zones.size() << "\n";
    std::cout << "Number of Entanglement Zones: " << config.entanglement_zones.size() << "\n";
    std::cout << "Number of AODs: " << config.aods.size() << "\n";

    // Print additional details as needed...

    return 0;
}
