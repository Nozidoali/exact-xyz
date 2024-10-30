#pragma once

#include "qstate.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace xyz
{

class QGate
{
private:
public:
  /* data */
  uint32_t target;
  QGate() = default;
  QGate( uint32_t target ) : target( target ){};
  virtual QState operator()( const QState& state, const bool reverse = false ) const = 0;
  virtual uint32_t num_cnots() const = 0;
  virtual std::string to_string() const = 0;
  friend std::ostream& operator<<( std::ostream& os, const QGate& obj );
};

class Rotation
{
public:
  /* data */
  static constexpr double eps = 1e-6;
  double theta;
  static bool is_trivial( double theta, bool use_x = false );
  Rotation( double theta ) : theta( theta ){};
};

class RY : public QGate, public Rotation
{
private:
public:
  /* data */
  RY( uint32_t target, double theta ) : Rotation( theta ), QGate( target ){};
  QState operator()( const QState& state, const bool reverse = false ) const;
  uint32_t num_cnots() const { return 0; };
  std::string to_string() const { return "ry(" + std::to_string( theta ) + ") q[" + std::to_string( target ) + "]"; };
};

class Controlled
{
private:
public:
  /* data */
  uint32_t ctrl;
  bool phase;
  Controlled( uint32_t ctrl, bool phase ) : ctrl( ctrl ), phase( phase ){};
};

class MultiControlled
{
private:
  /* data */
public:
  std::vector<uint32_t> ctrls;
  std::vector<bool> phases;
  MultiControlled( std::vector<uint32_t> ctrls ) : ctrls( ctrls ){};
};

class CRY : public Controlled, public RY
{
private:
public:
  /* data */
  CRY( uint32_t ctrl, bool phase, double theta, uint32_t target ) : Controlled( ctrl, phase ), RY( target, theta ){};
  QState operator()( const QState& state, const bool reverse = false ) const;
  std::string to_string() const { return "cry(" + std::to_string( theta ) + ") q[" + std::to_string( ctrl ) + "], q[" + std::to_string( target ) + "]"; };
  uint32_t num_cnots() const { return 2; };
};

class MCRY : public MultiControlled, public RY
{
private:
  /* data */
public:
  MCRY( std::vector<uint32_t> ctrls, double theta, uint32_t target ) : MultiControlled( ctrls ), RY( target, theta ){};
};

class X : public QGate
{
private:
  /* data */
public:
  using QGate::QGate;
  QState operator()( const QState& state, const bool reverse = false ) const;
  uint32_t num_cnots() const { return 0; };
  std::string to_string() const { return "x q[" + std::to_string( target ) + "]"; };
};

class CX : public Controlled, public X
{
private:
  /* data */
public:
  CX( uint32_t ctrl, bool phase, uint32_t target ) : Controlled( ctrl, phase ), X( target ){};
  QState operator()( const QState& state, const bool reverse = false ) const;
  std::string to_string() const { return "cx q[" + std::to_string( ctrl ) + "], q[" + std::to_string( target ) + "]"; };
  uint32_t num_cnots() const { return 1; };
};

} // namespace xyz