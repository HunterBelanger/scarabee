#ifndef SURFACE_H
#define SURFACE_H

#include <moc/vector.hpp>
#include <moc/direction.hpp>

#include <array>

class Surface {
 public:
  enum class Type : char { None, XPlane, YPlane, Plane, Cylinder };
  enum class Side : bool { Positive, Negative };

  Surface() : params_(), type_(Type::None) {}

  Side side(const Vector& r, const Direction& u) const;
  double distance(const Vector& r, const Direction& u) const;

  Type& type() { return type_; }
  const Type& type() const { return type_; }

  // Parameters for XPlane, YPlane, and Cylinder
  double& x0() { return params_[0]; }
  const double& x0() const { return params_[0]; }

  double& y0() { return params_[1]; }
  const double& y0() const { return params_[1]; }

  double& r() { return params_[2]; }
  const double& r() const { return params_[2]; }

  // Parameters for Plane
  double& A() { return params_[0]; }
  const double& A() const { return params_[0]; }

  double& B() { return params_[1]; }
  const double& B() const { return params_[1]; }

  double& C() { return params_[2]; }
  const double& C() const { return params_[2]; }

 private:
  std::array<double, 3> params_;
  Type type_;
};

#endif