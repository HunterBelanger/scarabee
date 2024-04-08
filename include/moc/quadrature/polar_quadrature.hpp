#ifndef POLAR_QUADRATURE_H
#define POLAR_QUADRATURE_H

#include <moc/quadrature/legendre.hpp>
#include <moc/quadrature/yamamoto_tabuchi.hpp>

#include <span>
#include <variant>

namespace scarabee {

using PolarQuadratureType =
    std::variant<Legendre<2>, Legendre<4>, Legendre<6>, Legendre<8>,
                 Legendre<10>, Legendre<12>, YamamotoTabuchi<2>,
                 YamamotoTabuchi<4>, YamamotoTabuchi<6>>;

class PolarQuadrature {
 public:
  template <typename P>
  PolarQuadrature(P pq) : pq_(pq), sin_(), invs_sin_(), wgt_() {
    sin_ = std::visit([](const auto& pq) { return pq.sin(); }, pq_);
    invs_sin_ = std::visit([](const auto& pq) { return pq.invs_sin(); }, pq_);
    wgt_ = std::visit([](const auto& pq) { return pq.wgt(); }, pq_);
  }

  const std::span<const double>& sin() const { return sin_; }
  const std::span<const double>& invs_sin() const { return invs_sin_; }
  const std::span<const double>& wgt() const { return wgt_; }

 private:
  PolarQuadratureType pq_;
  std::span<const double> sin_;
  std::span<const double> invs_sin_;
  std::span<const double> wgt_;
};

}  // namespace scarabee

#endif
