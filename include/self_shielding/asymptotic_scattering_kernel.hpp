#ifndef SCARABEE_ASYMPTOTIC_SCATTERING_KERNEL_H
#define SCARABEE_ASYMPTOTIC_SCATTERING_KERNEL_H

#include <utils/gauss_kronrod.hpp>
#include <utils/scarabee_exception.hpp>

#include <PapillonNDL/angle_distribution.hpp>

#include <cmath>
#include <utility>

/*
 * This is the asymptotic scattering kernel, which assumes that the target
 * nuclide is initially at rest.
 */
class AsymptoticScatteringKernel {
 public:
  AsymptoticScatteringKernel(const pndl::AngleDistribution& dist, double A,
                             double Q, double Elow, double Ehi, unsigned int l)
      : angle_dist_(dist), A_(A), Q_(Q), Eout_low_(Elow), Eout_hi_(Ehi), l_(l) {
    check_values();
  }

  double operator()(const double& E) const {
    // Since this is the asymptotic kernel, up scattering is impossible !
    // As such, we can imediately return 0 if E < Eout_low.
    if (E < Eout_low_) {
      return 0.;
    }

    // Get R constant for Center of Mass - Lab angle conversion
    const double R_ = R(E);

    // Get bounds of integration
    double omega_low = cm_cosine(E, Eout_low_, R_);
    double omega_hi = cm_cosine(E, Eout_hi_, R_);

    // If scattering angles are outside of [-1,1], we need to truncate the
    // domain of integration.
    if (omega_low < -1. && omega_hi < -1.) {
      return 0.;
    } else if (omega_low > 1. && omega_hi > 1.) {
      return 0.;
    } else if (omega_low < -1.) {
      omega_low = -1.;
    }
    if (omega_hi > 1.) {
      omega_hi = 1.;
    }

    // Create integrand function
    auto func = [this, omega_low, omega_hi, R_, E](const double& omega) {
      const double mu = lab_cosine(R_, omega);
      return this->angle_dist_.pdf(E, omega) * std::legendre(this->l_, mu);
    };

    // Integrate function
    GaussKronrodQuadrature<21> gk;
    auto I = gk.integrate(func, omega_low, omega_hi, 0.001, 100);

    // TODO check error on integral

    return I.first;
  }

  double A() const { return A_; }

  void set_A(const double& A) {
    A_ = A;
    check_values();
  }

  double Q() const { return Q_; }

  void set_Q(const double& Q) { Q_ = Q; }

  std::pair<double, double> Eout_bounds() const {
    return {Eout_low_, Eout_hi_};
  }

  void set_Eout_bounds(const std::pair<double, double>& Eout_bounds) {
    Eout_low_ = Eout_bounds.first;
    Eout_hi_ = Eout_bounds.second;
    check_values();
  }

  unsigned int legendre_order() const { return l_; }

  void set_legendre_order(const unsigned int& l) { l_ = l; }

 private:
  pndl::AngleDistribution
      angle_dist_;  // Scattering distribution in CM frame for reaction
  double A_, Q_;
  double Eout_low_, Eout_hi_;  // Outgoing energy group bounds
  unsigned int l_;             // Legendre order

  double R(const double& E) const {
    return A_ * std::sqrt(1. + (((A_ + 1.) * Q_) / (A_ * E)));
  }

  double lab_cosine(const double& R, const double& omega) const {
    return (1. + R * omega) / std::sqrt(1. + R * R + 2. * R * omega);
  }

  double cm_cosine(const double& E, const double& Eout, const double& R) const {
    const double num = (Eout * (1. + A_) * (1. + A_)) - (E * (1. + (R * R)));
    const double denom = 2. * R * E;
    return num / denom;
  }

  void check_values() const {
    if (A_ <= 0.) {
      throw ScarabeeException("A must be > 0.");
    }

    if (Eout_low_ <= 0.) {
      throw ScarabeeException("Eout_low must be > 0.");
    }

    if (Eout_hi_ <= 0.) {
      throw ScarabeeException("Eout_hi must be > 0.");
    }

    if (Eout_low_ >= Eout_hi_) {
      throw ScarabeeException("Eout_low must be < Eout_hi.");
    }
  }
};

#endif
