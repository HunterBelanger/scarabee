#ifndef YAMAMOTO_TABUCHI_H
#define YAMAMOTO_TABUCHI_H

#include <array>
#include <cmath>
#include <span>

namespace scarabee {

// Available N are 2, 4, and 6

// A. YAMAMOTO, M. TABUCHI, N. SUGIMURA, T. USHIO, and M. MORI, “Derivation of
// Optimum Polar Angle Quadrature Set for the Method of Characteristics Based on
// Approximation Error for the Bickley Function,” J. Nucl. Sci. Technol., vol.
// 44, no. 2, pp. 129–136, 2007, doi: 10.1080/18811248.2007.9711266.

template <std::size_t N>
class YamamotoTabuchi {
 public:
  std::span<const double> sin() const { return {sin_.begin(), sin_.end()}; }

  std::span<const double> invs_sin() const {
    return {invs_sin_.begin(), invs_sin_.end()};
  }

  std::span<const double> wgt() const { return {wgt_.begin(), wgt_.end()}; }

  std::span<const double> wsin() const { return {wsin_.begin(), wsin_.end()}; }

  std::span<const double> polar_angle() const {
    return {polar_angle_.begin(), polar_angle_.end()};
  }

 private:
  static inline const std::array<double, N / 2> invs_sin_;
  static inline const std::array<double, N / 2> wsin_;
  static inline const std::array<double, N / 2> sin_;
  static inline const std::array<double, N / 2> wgt_;
  static inline const std::array<double, N / 2> polar_angle_;
};

// A. YAMAMOTO, M. TABUCHI, N. SUGIMURA, T. USHIO, and M. MORI, “Derivation of
// Optimum Polar Angle Quadrature Set for the Method of Characteristics Based on
// Approximation Error for the Bickley Function,” J. Nucl. Sci. Technol., vol.
// 44, no. 2, pp. 129–136, 2007, doi: 10.1080/18811248.2007.9711266.

// N = 2
template <>
inline const std::array<double, 1> YamamotoTabuchi<2>::sin_ = {0.798184};
template <>
inline const std::array<double, 1> YamamotoTabuchi<2>::invs_sin_ = {1. /
                                                                    sin_[0]};
template <>
inline const std::array<double, 1> YamamotoTabuchi<2>::wgt_ = {1.000000};
template <>
inline const std::array<double, 1> YamamotoTabuchi<2>::wsin_ = {wgt_[0] *
                                                                sin_[0]};
template <>
inline const std::array<double, 1> YamamotoTabuchi<2>::polar_angle_ = {
    std::asin(sin_[0])};

// N = 4
template <>
inline const std::array<double, 2> YamamotoTabuchi<4>::sin_ = {0.363900,
                                                               0.899900};
template <>
inline const std::array<double, 2> YamamotoTabuchi<4>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1]};
template <>
inline const std::array<double, 2> YamamotoTabuchi<4>::wgt_ = {0.212854,
                                                               0.787146};
template <>
inline const std::array<double, 2> YamamotoTabuchi<4>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1]};
template <>
inline const std::array<double, 2> YamamotoTabuchi<4>::polar_angle_ = {
    std::asin(sin_[0]), std::asin(sin_[1])};

// N = 6
template <>
inline const std::array<double, 3> YamamotoTabuchi<6>::sin_ = {
    0.166648, 0.537707, 0.932954};
template <>
inline const std::array<double, 3> YamamotoTabuchi<6>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1], 1. / sin_[2]};
template <>
inline const std::array<double, 3> YamamotoTabuchi<6>::wgt_ = {
    0.046233, 0.283619, 0.670148};
template <>
inline const std::array<double, 3> YamamotoTabuchi<6>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1], wgt_[2] * sin_[2]};
template <>
inline const std::array<double, 3> YamamotoTabuchi<6>::polar_angle_ = {
    std::asin(sin_[0]), std::asin(sin_[1]), std::asin(sin_[2])};

}  // namespace scarabee

#endif
