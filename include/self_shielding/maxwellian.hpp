#ifndef SCARABEE_MAXWELLIAN_H
#define SCARABEE_MAXWELLIAN_H

#include <utils/scarabee_exception.hpp>

#include <cmath>

class Maxwellian {
  public:
    Maxwellian(double kT, double C = 1.): kT_{kT}, C_{C} {
      this->check_values(); 
    }

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
      this->check_values();
    }

    void set_temperature(double kT) {
      kT_ = kT;
      this->check_values();
    }

  private:
    double kT_; // Temperature in eV
    double C_;

    void check_values() const {
      if (kT_ <= 0.) {
        throw ScarabeeException("kT must be > 0.");
      }

      if (C_ <= 0.) {
        throw ScarabeeException("Constant must be > 0.");
      }
    }
};

#endif