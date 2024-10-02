#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

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
};

}  // namespace scarabee

#endif
