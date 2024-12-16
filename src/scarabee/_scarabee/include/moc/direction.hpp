#ifndef DIRECTION_H
#define DIRECTION_H

#include <moc/vector.hpp>
#include <utils/constants.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>

namespace scarabee {

//============================================================================
// Direction Class
//----------------------------------------------------------------------------
class Direction : public Vector {
 public:
  Direction() : Vector(1., 0.) {}
  Direction(double i_x, double i_y) : Vector(i_x, i_y) {
    double magnitude = this->norm();
    this->x_ /= magnitude;
    this->y_ /= magnitude;
  }
  Direction(double phi) : Vector(0., 1.) {
    if (phi < 0.)
      phi = 0.;
    else if (phi > 2 * PI)
      phi = 2 * PI;

    this->x_ = std::cos(phi);
    this->y_ = std::sin(phi);
  }
  ~Direction() = default;

  Direction operator-() const {
    Direction d;
    d.x_ = -this->x_;
    d.y_ = -this->y_;
    return d;
  }

 private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<Vector>(this));
  }
};  // Direction

//============================================================================
// Addition Operators
inline Vector operator+(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() + d2.x(), d1.y() + d2.y());
}

inline Vector operator+(const Direction& d, const Vector& v) {
  return Vector(d.x() + v.x(), d.y() + v.y());
}

inline Vector operator+(const Vector& v, const Direction& d) {
  return Vector(d.x() + v.x(), d.y() + v.y());
}

//============================================================================
// Subtraction Operators
inline Vector operator-(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() - d2.x(), d1.y() - d2.y());
}

inline Vector operator-(const Direction& d, const Vector& v) {
  return Vector(d.x() - v.x(), d.y() - v.y());
}

inline Vector operator-(const Vector& v, const Direction& d) {
  return Vector(v.x() - d.x(), v.y() - d.y());
}

//============================================================================
// Dot Product Operators
inline double operator*(const Direction& d, const Direction& v) {
  return d.dot(v);
}

inline double operator*(const Direction& d, const Vector& v) {
  return d.dot(v);
}

inline double operator*(const Vector& v, const Direction& d) {
  return d.dot(v);
}

//============================================================================
// Scaling Operators
inline Vector operator*(const Direction& d, double c) {
  return Vector(d.x() * c, d.y() * c);
}

inline Vector operator*(double c, const Direction& d) {
  return Vector(d.x() * c, d.y() * c);
}

inline Vector operator/(const Direction& d, double c) {
  return Vector(d.x() / c, d.y() / c);
}

inline std::ostream& operator<<(std::ostream& output, const Direction& d) {
  output << "<<" << d.x() << "," << d.y() << ">>";
  return output;
}

}  // namespace scarabee

#endif
