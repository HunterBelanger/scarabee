#ifndef SCARABEE_WATT_H
#define SCARABEE_WATT_H

#include <cmath>

class Watt {
  public:
    Watt(double a, double b, double C = 1.): a_{a}, b_{b}, C_{C} {}

    double operator()(double E) const {
      return C_ * std::exp(-E/a_) * std::sinh(std::sqrt(b_ * E));
    }

    double a() const {
      return a_;
    }

    double b() const {
      return b_;
    }

    double constant() const {
      return C_;
    }

    void set_a(double a) {
      a_ = a;
    }

    void set_b(double b) {
      b_ = b;
    }

    void set_constant(double C) {
      C_ = C;
    }

  private:
    double a_;
    double b_;
    double C_;
};

#endif