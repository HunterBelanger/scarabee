#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

//============================================================================
// Vector Class
//----------------------------------------------------------------------------
class Vector {
 public:
  Vector(double i_x, double i_y) : x_(i_x), y_(i_y) {};
  ~Vector() = default;

  double x() const { return x_; }
  double y() const { return y_; }

  double dot(const Vector& v) const {
    return x_ * v.x() + y_ * v.y();
  }

  double norm() const { return std::sqrt(x_ * x_ + y_ * y_); }

 protected:
  double x_, y_;
};

//============================================================================
// Overloaded Operator Declarations
//----------------------------------------------------------------------------
inline Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.x() + v2.x(), v1.y() + v2.y());
}

inline Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.x() - v2.x(), v1.y() - v2.y());
}

inline Vector operator*(const Vector& v, double d) {
  return Vector(v.x() * d, v.y() * d);
}

inline Vector operator*(double d, const Vector& v) { return v * d; }

inline Vector operator/(const Vector& v, double d) {
  return Vector(v.x() / d, v.y() / d);
}

inline double operator*(const Vector& v1, const Vector& v2) {
  return v1.dot(v2);
}

inline std::ostream& operator<<(std::ostream& output, const Vector& v) {
  output << "<" << v.x() << "," << v.y() << ">";
  return output;
}

#endif
