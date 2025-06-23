#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <utils/serialization.hpp>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

#include <cereal/cereal.hpp>

#include <cstdint>
#include <span>

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

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(all_harmonics_), CEREAL_NVP(L_), CEREAL_NVP(Nlj_));
  }
};

}  // namespace scarabee

#endif
