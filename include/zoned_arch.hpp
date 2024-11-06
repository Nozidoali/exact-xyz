#pragma once

#include "arch.hpp"
#include "qcircuit.hpp"

namespace xyz
{

namespace zoned
{
class Gate {
public:
    Gate() = default;
};

class Gate1Q : public Gate {
public:
    std::string name; // type of the U3 gate
    int q;
    Gate1Q(const std::string& name, int q) : name(name), q(q) {}
};

class Gate2Q : public Gate {
public:
    int id; // Gate ID
    int q0, q1;
    Gate2Q(int id, int q0, int q1) : id(id), q0(q0), q1(q1) {}
};

class Inst {
public:
  int id;
  std::string name;
  double begin_time, end_time;

  Inst(int id, const std::string& name, double begin, double end)
    : id(id), name(name), begin_time(begin), end_time(end) {}

  virtual ~Inst() = default;
  virtual inline double duration() const = 0; // Virtual method to return the duration of the instruction
  virtual void print() const = 0; // Virtual method to print details of the instruction
};

// InitInst class
class InitInst : public Inst {
public:
  std::vector<std::vector<int>> init_locs;

  InitInst(int id, double begin, double end, std::vector<std::vector<int>> locs)
    : Inst(id, "InitInst", begin, end), init_locs(locs) {}

  void print() const override {
    std::cout << name << ": id=" << id << ", begin_time=" << begin_time
          << ", end_time=" << end_time << ", init_locs.size=" << init_locs.size() << std::endl;
  }

  inline double duration() const override {
    return 0.0;
  }
};

// OneQGateInst class
class OneQGateInst : public Inst {
public:
  std::string unitary;
  std::vector<std::vector<int>> locs;
  std::vector<std::pair<std::string, int>> gates;
  std::map<std::string, std::vector<int>> dependency;

  OneQGateInst(int id, std::string unitary, double begin, double end,
         std::vector<std::vector<int>> locs,
         std::vector<std::pair<std::string, int>> gates,
         std::map<std::string, std::vector<int>> dependency)
    : Inst(id, "OneQGateInst", begin, end), unitary(unitary), 
      locs(locs), gates(gates), dependency(dependency) {}

  void print() const override {
    std::cout << name << ": id=" << id << ", unitary=" << unitary
          << ", begin_time=" << begin_time << ", end_time=" << end_time
          << ", locs.size=" << locs.size() << ", gates.size=" << gates.size() << std::endl;
  }

  inline double duration() const override {
    return constants::t_gate1q;
  }
};

// RearrangeJobInst class
class RearrangeJobInst : public Inst {
public:
  int aod_id;
  std::vector<int> aod_qubits;
  std::vector<std::vector<int>> begin_locs, end_locs;
  std::map<std::string, std::vector<int>> dependency;

  RearrangeJobInst(int id, int aod_id, double begin, double end,
           std::vector<int> aod_qubits, std::vector<std::vector<int>> begin_locs,
           std::vector<std::vector<int>> end_locs,
           std::map<std::string, std::vector<int>> dependency)
    : Inst(id, "RearrangeJobInst", begin, end), aod_id(aod_id), 
      aod_qubits(aod_qubits), begin_locs(begin_locs), end_locs(end_locs), dependency(dependency) {}

  void print() const override {
    std::cout << name << ": id=" << id << ", aod_id=" << aod_id
          << ", begin_time=" << begin_time << ", end_time=" << end_time
          << ", aod_qubits.size=" << aod_qubits.size() << std::endl;
  }

  inline double duration() const override {
    return constants::t_transfer;
  }
};

// RydbergInst class
class RydbergInst : public Inst {
public:
  int zone_id;
  std::vector<std::shared_ptr<Gate>> gates;
  std::map<std::string, std::vector<int>> dependency;

  RydbergInst(int id, int zone_id, double begin, double end,
        std::vector<std::shared_ptr<Gate>> gates,
        std::map<std::string, std::vector<int>> dependency)
    : Inst(id, "RydbergInst", begin, end), zone_id(zone_id), 
      gates(gates), dependency(dependency) {}

  void print() const override {
    std::cout << name << ": id=" << id << ", zone_id=" << zone_id
          << ", begin_time=" << begin_time << ", end_time=" << end_time
          << ", gates.size=" << gates.size() << std::endl;
  }

  inline double duration() const override {
    return constants::t_rydberg;
  }
};

using Schedule = std::vector<std::shared_ptr<Inst>>;

Schedule layout_synthesis( const QCircuit& circuit, const Config& config );
Schedule load_schedule( const std::string& filename );
void simulate( const Schedule& schedule, const Config& config );
} // namespace zoned
} // namespace xyz