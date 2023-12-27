#ifndef SCARABEE_WATT_H
#define SCARABEE_WATT_H

#include <cmath>

class Watt {
 public:
  Watt(double a, double b, double C = 1.) : a_{a}, b_{b}, C_{C} {
    this->check_values();
  }

  double operator()(double E) const {
    return C_ * std::exp(-E / a_) * std::sinh(std::sqrt(b_ * E));
  }

  double a() const { return a_; }

  double b() const { return b_; }

  double constant() const { return C_; }

  void set_a(double a) {
    a_ = a;
    this->check_values();
  }

  void set_b(double b) {
    b_ = b;
    this->check_values();
  }

  void set_constant(double C) {
    C_ = C;
    this->check_values();
  }

 private:
  double a_;
  double b_;
  double C_;

  void check_values() const {
    if (a_ <= 0.) {
      throw ScarabeeException("a must be > 0.");
    }

    if (b_ <= 0.) {
      throw ScarabeeException("b must be > 0.");
    }

    if (C_ <= 0.) {
      throw ScarabeeException("Constant must be > 0.");
    }
  }
};

#endif