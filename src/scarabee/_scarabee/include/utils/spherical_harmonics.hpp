#ifndef SPHERICAL_HARMONICS_V2_H
#define SPHERICAL_HARMONICS_V2_H

#include <utils/constants.hpp>

#include <cstdint>
#include <span>

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace scarabee {

class SphericalHarmonics {
 public:
  SphericalHarmonics() = default;

  SphericalHarmonics(const std::size_t& L,
                     const std::vector<double>& azimuthal_angle,
                     const std::vector<double>& polar_angle);

  std::span<const double> spherical_harmonics(
      const std::size_t phi_index, const std::size_t theta_index) const {
    std::span<const double> Ylj(&all_harmonics_(phi_index, theta_index, 0),
                                Nlj_);
    return Ylj;
  }

 private:
  xt::xtensor<double, 3> all_harmonics_;  // Phi (azimuthal angle), Theta (polar
                                          // angle), l/j spherical harmonics
  std::size_t L_;
  std::size_t Nlj_;

  // general method to evaluate the spherical harmonics
  double eval_spherical_harmonics(const std::size_t& l, const int& j,
                                  const double& theta, const double& phi) const;

  // factorial
  double factorial(const std::size_t& N) const {
    if (N == 1 || N == 0)
      return 1.;
    else
      return static_cast<double>(N) * factorial(N - 1);
  }

  // method of evaluae the legendre
  double eval_legendre(const std::size_t& order, const double& x) const {
    if (order > 1) {
      return (static_cast<double>(2 * order - 1) * eval_legendre(order - 1, x) *
                  x -
              static_cast<double>(order - 1) * eval_legendre(order - 2, x)) /
             static_cast<double>(order);
    } else if (order == 1) {
      return x;
    }

    // if order == 0
    return 1.;
  }

  // method to get differentiation of legendre-polynomial
  double eval_legendre_differentiation(const std::size_t& order,
                                       const std::size_t& n,
                                       const double& x) const;
  // method to get the associated legendre
  double eval_associated_legendre(const std::size_t& order, const int& j,
                                  const double& x) const;
};

}  // namespace scarabee

#endif
