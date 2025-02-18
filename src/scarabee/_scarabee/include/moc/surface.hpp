#ifndef SURFACE_H
#define SURFACE_H

#include <moc/vector.hpp>
#include <moc/direction.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/base_class.hpp>

#include <array>

namespace scarabee {

class Surface {
 public:
  enum class Type : char { XPlane, YPlane, Plane, Cylinder };
  enum class Side : bool { Positive, Negative };

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

 protected:
  Surface(Type type) : params_(), type_(type) {}

 private:
  Surface() : params_(), type_(Type::XPlane) {}

  std::array<double, 3> params_;
  Type type_;

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(params_), CEREAL_NVP(type_));
  }
};

class XPlane : public Surface {
 public:
  XPlane(double x0) : Surface(Type::XPlane) { this->x0() = x0; }

 private:
  XPlane(): Surface(Type::XPlane) {}

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<Surface>(this));
  }
};

class YPlane : public Surface {
 public:
  YPlane(double y0) : Surface(Type::YPlane) { this->y0() = y0; }

 private:
  YPlane(): Surface(Type::YPlane) {}

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<Surface>(this));
  }
};

class Plane : public Surface {
 public:
  Plane(double A, double B, double C) : Surface(Type::Plane) {
    this->A() = A;
    this->B() = B;
    this->C() = C;
  }

 private:
  Plane(): Surface(Type::Plane) {}

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<Surface>(this));
  }
};

class Cylinder : public Surface {
 public:
  Cylinder(double x0, double y0, double r) : Surface(Type::Cylinder) {
    this->x0() = x0;
    this->y0() = y0;
    this->r() = r;
  }

 private:
  Cylinder(): Surface(Type::Cylinder) {}

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<Surface>(this));
  }
};

}  // namespace scarabee

#endif
