#ifndef SCARABEE_TRANSPORT_CROSS_SECTIONS_H
#define SCARABEE_TRANSPORT_CROSS_SECTIONS_H

#include <utils/scarabee_exception.hpp>

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <cstdint>

namespace scarabee {

class TransportXS {
 public:
  TransportXS();

  std::size_t ngroups() const { return Et_.size(); }

  bool fissile() const { return fissile_; }

  const xt::xtensor<double, 1>& Et() const { return Et_; }

  double Et(std::size_t g) const { return Et_(g); }

  double Ea(std::size_t g) const { return Ea_(g); }

  double Ef(std::size_t g) const {
    if (fissile_) return Ef_(g);
    return 0.;
  }

  double Er(std::size_t g) const { return Ea(g) + Es(g) - Es(g, g); }

  double nu(std::size_t g) const {
    if (fissile_) return nu_(g);
    return 0.;
  }

  double chi(std::size_t g) const {
    if (fissile_) return chi_(g);
    return 0.;
  }

  double Es(std::size_t gin, std::size_t gout) const { return Es_(gin, gout); }

  double Es(std::size_t gin) const {
    return xt::sum(xt::view(Es_, gin, xt::all()))();
  }

  // Operators for constructing compound cross sections
  TransportXS operator+(const TransportXS& R) const;
  TransportXS operator*(double N) const;
  TransportXS& operator+=(const TransportXS& R);
  TransportXS& operator*=(double N);

  // private:
  xt::xtensor<double, 2> Es_;   // Scattering matrix
  xt::xtensor<double, 1> Et_;   // Total xs
  xt::xtensor<double, 1> Ea_;   // Absorption xs
  xt::xtensor<double, 1> Ef_;   // Fission xs
  xt::xtensor<double, 1> nu_;   // Fission yields
  xt::xtensor<double, 1> chi_;  // Fission spectrum
  bool fissile_;                // Fissile indicator
};

}  // namespace scarabee

#endif
