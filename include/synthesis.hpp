#pragma once

#include "arch.hpp"
#include "qcircuit.hpp"

namespace xyz
{
class Inst {
public:
    int id;
    std::string name;

    Inst(int id, const std::string& name) : id(id), name(name) {}

    virtual ~Inst() = default;
    virtual void print() const = 0; // Virtual method to print details of the instruction
};

// InitInst class
class InitInst : public Inst {
public:
    double begin_time, end_time;
    std::vector<std::vector<int>> init_locs;

    InitInst(int id, double begin, double end, std::vector<std::vector<int>> locs)
        : Inst(id, "InitInst"), begin_time(begin), end_time(end), init_locs(locs) {}

    void print() const override {
        std::cout << name << ": id=" << id << ", begin_time=" << begin_time
                  << ", end_time=" << end_time << ", init_locs.size=" << init_locs.size() << std::endl;
    }
};

// OneQGateInst class
class OneQGateInst : public Inst {
public:
    std::string unitary;
    double begin_time, end_time;
    std::vector<std::vector<int>> locs;
    std::vector<std::pair<std::string, int>> gates;
    std::map<std::string, std::vector<int>> dependency;

    OneQGateInst(int id, std::string unitary, double begin, double end,
                 std::vector<std::vector<int>> locs,
                 std::vector<std::pair<std::string, int>> gates,
                 std::map<std::string, std::vector<int>> dependency)
        : Inst(id, "OneQGateInst"), unitary(unitary), begin_time(begin), end_time(end), 
          locs(locs), gates(gates), dependency(dependency) {}

    void print() const override {
        std::cout << name << ": id=" << id << ", unitary=" << unitary
                  << ", begin_time=" << begin_time << ", end_time=" << end_time
                  << ", locs.size=" << locs.size() << ", gates.size=" << gates.size() << std::endl;
    }
};

// RearrangeJobInst class
class RearrangeJobInst : public Inst {
public:
    int aod_id;
    double begin_time, end_time;
    std::vector<int> aod_qubits;
    std::vector<std::vector<int>> begin_locs, end_locs;
    std::map<std::string, std::vector<int>> dependency;

    RearrangeJobInst(int id, int aod_id, double begin, double end,
                     std::vector<int> aod_qubits, std::vector<std::vector<int>> begin_locs,
                     std::vector<std::vector<int>> end_locs,
                     std::map<std::string, std::vector<int>> dependency)
        : Inst(id, "RearrangeJobInst"), aod_id(aod_id), begin_time(begin), end_time(end), 
          aod_qubits(aod_qubits), begin_locs(begin_locs), end_locs(end_locs), dependency(dependency) {}

    void print() const override {
        std::cout << name << ": id=" << id << ", aod_id=" << aod_id
                  << ", begin_time=" << begin_time << ", end_time=" << end_time
                  << ", aod_qubits.size=" << aod_qubits.size() << std::endl;
    }
};

// RydbergInst class
class RydbergInst : public Inst {
public:
    int zone_id;
    double begin_time, end_time;
    std::vector<std::pair<int, int>> gates;
    std::map<std::string, std::vector<int>> dependency;

    RydbergInst(int id, int zone_id, double begin, double end,
                std::vector<std::pair<int, int>> gates,
                std::map<std::string, std::vector<int>> dependency)
        : Inst(id, "RydbergInst"), zone_id(zone_id), begin_time(begin), end_time(end), 
          gates(gates), dependency(dependency) {}

    void print() const override {
        std::cout << name << ": id=" << id << ", zone_id=" << zone_id
                  << ", begin_time=" << begin_time << ", end_time=" << end_time
                  << ", gates.size=" << gates.size() << std::endl;
    }
};

using Schedule = std::vector<std::shared_ptr<Inst>>;

Schedule layout_synthesis( const QCircuit& circuit, const Config& config );

Schedule load_schedule( const std::string& filename );
void simulate( const Schedule& schedule );

} // namespace xyz