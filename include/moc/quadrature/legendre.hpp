#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <array>
#include <span>

namespace scarabee {

// Available N are 2, 4, 6, 8, 10, and 12

template <std::size_t N>
class Legendre {
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
