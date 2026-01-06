#pragma once

#include "qstate.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

namespace xyz {
using namespace std::complex_literals;

namespace constants {
constexpr double sqrt2     = 1.4142135623730951;
constexpr double sqrt2_inv = 0.7071067811865475;
} // namespace constants

class QGate {
  public:
    uint32_t target;
    QGate() = default;
    QGate(uint32_t target) : target(target) {};
    virtual QRState     operator()(const QRState& state, const bool reverse = false) const = 0;
    virtual uint32_t    num_cnots() const                                                  = 0;
    virtual std::string to_string() const                                                  = 0;

    friend std::ostream& operator<<(std::ostream& os, const QGate& obj) {
        os << obj.to_string();
        return os;
    };
    virtual std::vector<uint32_t> qbits() const { return {target}; };
};

class Rotation {
  public:
    static constexpr double eps = 1e-6;
    double                  theta;
    static bool             is_trivial(double theta, bool use_x = false);
    Rotation(double theta) : theta(theta) {};
};

class RU2 {
  public:
    std::array<double, 2> c00;
    std::array<double, 2> c01;
    std::array<double, 2> c10;
    std::array<double, 2> c11;

    RU2(double _c00, double _c01, double _c10, double _c11) {
        c00 = {_c00, _c00};
        c01 = {_c01, _c10};
        c10 = {_c10, _c01};
        c11 = {_c11, _c11};
    };
    QRState operator()(const QRState& state, uint32_t target, const bool reverse = false) const;
};

class U2 {
  public:
    std::array<std::complex<double>, 2> c00;
    std::array<std::complex<double>, 2> c01;
    std::array<std::complex<double>, 2> c10;
    std::array<std::complex<double>, 2> c11;

    U2(std::complex<double> _c00, std::complex<double> _c01, std::complex<double> _c10, std::complex<double> _c11) {
        c00 = {_c00, std::conj(_c00)};
        c01 = {_c01, std::conj(_c10)};
        c10 = {_c10, std::conj(_c01)};
        c11 = {_c11, std::conj(_c11)};
    };
    QRState operator()(const QRState& state, uint32_t target, const bool reverse = false) const;
    QState  operator()(const QState& state, uint32_t target, const bool reverse = false) const;
};

class RY : public QGate, public Rotation, public RU2 {
  public:
    RY(uint32_t target, double theta)
        : Rotation(theta), QGate(target), RU2(cos(theta / 2), sin(theta / 2), -sin(theta / 2), cos(theta / 2)) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        return RU2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override {
        return "ry(" + std::to_string(theta) + ") q[" + std::to_string(target) + "]";
    };
};

class X : public QGate {
  public:
    using QGate::QGate;
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "x q[" + std::to_string(target) + "]"; };
};

class H : public QGate, public RU2 {
  public:
    using QGate::QGate;
    H(uint32_t target)
        : QGate(target),
          RU2(constants::sqrt2_inv, constants::sqrt2_inv, constants::sqrt2_inv, -constants::sqrt2_inv) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        return RU2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "h q[" + std::to_string(target) + "]"; };
};

class T : public QGate, public U2 {
  public:
    using QGate::QGate;
    T(uint32_t target) : QGate(target), U2(1, 0, 0, std::exp(1i * M_PI / 4.0)) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        return U2::operator()(state, target, reverse);
    };
    QState operator()(const QState& state, const bool reverse = false) const {
        return U2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "t q[" + std::to_string(target) + "]"; };
};

class Tdg : public QGate, public U2 {
  public:
    using QGate::QGate;
    Tdg(uint32_t target) : QGate(target), U2(1, 0, 0, std::exp(-1i * M_PI / 4.0)) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        return U2::operator()(state, target, reverse);
    };
    QState operator()(const QState& state, const bool reverse = false) const {
        return U2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "tdg q[" + std::to_string(target) + "]"; };
};

class Sdg : public QGate, public U2 {
  public:
    using QGate::QGate;
    Sdg(uint32_t target) : QGate(target), U2(1, 0, 0, std::exp(-1i * M_PI / 2.0)) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        throw std::runtime_error("Sdg introduces complex phase; use QState simulation");
    };
    QState operator()(const QState& state, const bool reverse = false) const {
        return U2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "sdg q[" + std::to_string(target) + "]"; };
};

class Z : public QGate, public RU2 {
  public:
    using QGate::QGate;
    Z(uint32_t target) : QGate(target), RU2(1, 0, 0, -1) {};
    QRState operator()(const QRState& state, const bool reverse = false) const override {
        return RU2::operator()(state, target, reverse);
    };
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "z q[" + std::to_string(target) + "]"; };
};

class S : public QGate, public U2 {
  public:
    using QGate::QGate;
    S(uint32_t target) : QGate(target), U2(1, 0, 0, 1i) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    uint32_t    num_cnots() const override { return 0; };
    std::string to_string() const override { return "s q[" + std::to_string(target) + "]"; };
};

class Controlled {
  public:
    uint32_t ctrl;
    bool     phase;
    Controlled(uint32_t ctrl, bool phase) : ctrl(ctrl), phase(phase) {};
    virtual std::vector<uint32_t> qbits() const { return {ctrl}; };
};

class MultiControlled {
  public:
    std::vector<uint32_t> ctrls;
    std::vector<bool>     phases;
    MultiControlled() = default;
    MultiControlled(std::vector<uint32_t> ctrls) : ctrls(ctrls), phases(ctrls.size(), true) {};
    MultiControlled(std::vector<uint32_t> ctrls, std::vector<bool> phases) : ctrls(ctrls), phases(phases) {};
    std::vector<uint32_t> qbits() const { return ctrls; };
};

class CRY : public Controlled, public RY {
  public:
    CRY(uint32_t ctrl, bool phase, double theta, uint32_t target) : Controlled(ctrl, phase), RY(target, theta) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    std::string to_string() const override {
        std::string gate = phase ? "cry" : "cry_false";
        return gate + "(" + std::to_string(theta) + ") q[" + std::to_string(ctrl) + "], q[" + std::to_string(target) +
               "]";
    };
    uint32_t              num_cnots() const override { return 2; };
    std::vector<uint32_t> qbits() const override { return {ctrl, target}; };
};

class MCRY : public MultiControlled, public RY {
  public:
    MCRY(std::vector<uint32_t> ctrls, double theta, uint32_t target) : MultiControlled(ctrls), RY(target, theta) {};
    MCRY(std::vector<uint32_t> ctrls, std::vector<bool> phases, double theta, uint32_t target)
        : MultiControlled(ctrls, phases), RY(target, theta) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    std::string to_string() const override {
        std::string gate = "mcry";
        for (uint32_t i = 0; i < ctrls.size(); i++)
            gate += "[" + std::to_string(ctrls[i]) + "]";
        return gate + "(" + std::to_string(theta) + ") q[" + std::to_string(target) + "]";
    };
    uint32_t num_cnots() const override { return 1 << ctrls.size(); };
};

class CX : public Controlled, public X {
  public:
    CX(uint32_t ctrl, bool phase, uint32_t target) : Controlled(ctrl, phase), X(target) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    std::string to_string() const override {
        std::string gate = phase ? "cx" : "cx_false";
        return gate + " q[" + std::to_string(ctrl) + "], q[" + std::to_string(target) + "]";
    };
    uint32_t              num_cnots() const override { return 1; };
    std::vector<uint32_t> qbits() const override { return {ctrl, target}; };
};

class CCX : public MultiControlled, public X {
  public:
    CCX(uint32_t ctrl1, uint32_t ctrl2, uint32_t target) : MultiControlled({ctrl1, ctrl2}), X(target) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    std::string to_string() const override {
        return "ccx q[" + std::to_string(ctrls[0]) + "], q[" + std::to_string(ctrls[1]) + "], q[" +
               std::to_string(target) + "]";
    };
    uint32_t              num_cnots() const override { return 2; };
    std::vector<uint32_t> qbits() const override { return {ctrls[0], ctrls[1], target}; };
};

class QROM_MCRY : public MultiControlled, public RY {
  public:
    double              eps;
    std::vector<double> rotation_table;
    QROM_MCRY(std::vector<uint32_t> ctrls, std::vector<bool> phases, std::vector<double> rotation_table,
              uint32_t target, double eps)
        : MultiControlled(ctrls, phases), RY(target, 0.0), eps(eps), rotation_table(rotation_table) {};
    QRState     operator()(const QRState& state, const bool reverse = false) const override;
    std::string to_string() const override {
        return "qrom_mcry q[" + std::to_string(target) + "], eps=" + std::to_string(eps);
    };
    uint32_t num_cnots() const override {
        uint32_t k = ctrls.size();
        return std::min(1u << k, (uint32_t)(k * k * std::log(1.0 / eps)));
    };
};

std::vector<std::shared_ptr<QGate>> decompose_mcry(const MCRY& gate);

} // namespace xyz
