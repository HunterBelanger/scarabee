#ifndef POLAR_QUADRATURE_H
#define POLAR_QUADRATURE_H

#include <moc/quadrature/legendre.hpp>
#include <moc/quadrature/yamamoto_tabuchi.hpp>

#include <span>
#include <variant>

using PolarQuadratureType =
    std::variant<Legendre<2>, Legendre<4>, Legendre<6>, Legendre<8>,
                 Legendre<10>, Legendre<12>, YamamotoTabuchi<2>,
                 YamamotoTabuchi<4>, YamamotoTabuchi<6>>;

class PolarQuadrature {
 public:
  template <typename P>
  PolarQuadrature(P pq) : pq_(pq), abscissae_(), weights_() {
    abscissae_ = std::visit([](const auto& pq) { return pq.abscissae(); }, pq_);
    weights_ = std::visit([](const auto& pq) { return pq.weights(); }, pq_);
  }

  const std::span<const double>& abscissae() const { return abscissae_; }
  const std::span<const double>& weights() const { return weights_; }

 private:
  PolarQuadratureType pq_;
  std::span<const double> abscissae_;
  std::span<const double> weights_;
};

#endif