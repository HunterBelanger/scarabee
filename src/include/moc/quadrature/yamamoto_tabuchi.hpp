#ifndef YAMAMOTO_TABUCHI_H
#define YAMAMOTO_TABUCHI_H

#include <array>
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

 private:
  static const std::array<double, N / 2> sin_;
  static const std::array<double, N / 2> invs_sin_;
  static const std::array<double, N / 2> wgt_;
};

}  // namespace scarabee

#endif
