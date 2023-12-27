#ifndef SCARABEE_ONE_OVER_E_H
#define SCARABEE_ONE_OVER_E_H

#include <utils/scarabee_exception.hpp>

class OneOverE {
 public:
  OneOverE(double C = 1.) : C_{C} { this->check_values(); }

  double operator()(double E) const { return C_ / E; }

  double constant() const { return C_; }

  void set_constant(double C) {
    C_ = C;
    this->check_values();
  }

 private:
  double C_;

  void check_values() const {
    if (C_ <= 0.) {
      throw ScarabeeException("Constant must be > 0.");
    }
  }
};

#endif