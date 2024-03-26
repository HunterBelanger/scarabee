#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <array>
#include <cstdint>
#include <span>

// Available N are 2, 4, 6, 8, 10, and 12

template <std::size_t N>
class Legendre {
 public:
  std::span<const double> abscissae() const {
    return {abscissae_.begin(), abscissae_.end()};
  }
  std::span<const double> weights() const {
    return {weights_.begin(), weights_.end()};
  }

 private:
  static const std::array<double, N / 2> abscissae_;
  static const std::array<double, N / 2> weights_;
};

#endif