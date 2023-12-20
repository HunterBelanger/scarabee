#ifndef SCARABEE_MAXWELLIAN_H
#define SCARABEE_MAXWELLIAN_H

#include <cmath>

class Maxwellian {
  public:
    Maxwellian(double kT, double C = 1.): kT_{kT}, C_{C} {}

    double operator()(double E) const {
      return C_ * std::sqrt(E) * std::exp(-E / kT_);
    }

    double constant() const {
      return C_;
    }

    double temperature() const {
      return kT_;
    }

    void set_constant(double C) {
      C_ = C;
    }

    void set_temperature(double kT) {
      kT_ = kT;
    }

  private:
    double kT_; // Temperature in eV
    double C_;
};

#endif