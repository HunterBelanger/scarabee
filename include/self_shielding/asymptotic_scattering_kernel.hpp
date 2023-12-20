#ifndef SCARABEE_ASYMPTOTIC_SCATTERING_KERNEL_H
#define SCARABEE_ASYMPTOTIC_SCATTERING_KERNEL_H

#include <utils/gauss_kronrod.hpp>

#include <PapillonNDL/angle_distribution.hpp>

#include <cmath>
#include <exception>

double R(const double& E, const double& A, const double& Q) {
  return A * std::sqrt(1. + (((A+1.)*Q)/(A*E)));
}

double lab_cosine(const double& R, const double& omega) {
  return (1. + R*omega) / std::sqrt(1. + R*R + 2.*R*omega);
}

double cm_cosine(const double& E, const double& Eout, const double& A, const double& R) {
  const double num = (Eout * (1. + A) * (1. + A)) - (E * (1. + (R * R)));
  const double denom = 2. * R * E;
  return num / denom;
}

/*
 * This is the asymptotic scattering kernel, which assumes that the target
 * nuclide is initially at rest.
 */
struct AsymptoticScatteringKernel {

  double operator()(const double& E) const {
    // It is required that the outgoing energy group bounds be in order
    if (Eout_low >= Eout_hi) {
        throw std::runtime_error("Eout_low must be < Eout_hi.");
    }

    // Since this is the asymptotic kernel, up scattering is impossible !
    // As such, we can imediately return 0 if E < Eout_low.
    if (E < Eout_low) { return 0.; }

    // Get R constant for Center of Mass - Lab angle conversion
    const double R_ = R(E, A, Q);

    // Get bounds of integration
    double omega_low = cm_cosine(E, Eout_low, A, R_);
    double omega_hi = cm_cosine(E, Eout_hi, A, R_);

    // If scattering angles are outside of [-1,1], we need to truncate the
    // domain of integration.
    if (omega_low < -1. && omega_hi < -1.) { return 0.; }
    else if (omega_low > 1. && omega_hi > 1.) { return 0.; }
    else if (omega_low < -1.) { omega_low = -1.; }
    if (omega_hi > 1.) { omega_hi = 1.; }

    // Create integrand function
    auto func = [this, omega_low, omega_hi, R_, E](const double& omega) {
      const double mu = lab_cosine(R_, omega);
      return this->angle_dist_.pdf(E, omega) * std::legendre(this->l, mu);
    };

    // Integrate function
    GaussKronrodQuadrature<21> gk;
    auto I = gk.integrate(func, omega_low, omega_hi, 0.001, 100);

    // TODO check error on integral

    return I.first;
  }

  pndl::AngleDistribution angle_dist_; // Scattering distribution in CM frame for reaction
  double A, Q;
  double Eout_low, Eout_hi; // Outgoing energy group bounds
  unsigned int l; // Legendre order
};

#endif