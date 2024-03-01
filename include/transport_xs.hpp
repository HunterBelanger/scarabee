#ifndef SCARABEE_TRANSPORT_CROSS_SECTIONS_H
#define SCARABEE_TRANSPORT_CROSS_SECTIONS_H

#include <utils/scarabee_exception.hpp>

#include <xtensor/xarray.hpp>

#include <cstdint>
#include <vector>

class TransportXS {
 public:
  TransportXS();

  std::size_t ngroups() const { return Et_.size(); }

  bool fissile() const { return fissile_; }

  double Et(std::uint32_t g) const { return Et_[g]; }

  double Ea(std::uint32_t g) const { return Ea_[g]; }

  double Ef(std::uint32_t g) const {
    if (fissile_) return Ef_[g];
    return 0.;
  }

  double Er(std::uint32_t g) const { return Ea(g) + Es(g) - Es(g, g); }

  double nu(std::uint32_t g) const {
    if (fissile_) return nu_[g];
    return 0.;
  }

  double chi(std::uint32_t g) const {
    if (fissile_) return chi_[g];
    return 0.;
  }

  double Es(std::uint32_t gin, std::uint32_t gout) const {
    return Es_(gin, gout);
  }

  double Es(std::uint32_t gin) const {
    double Es_ret = 0.;
    for (std::uint32_t gout = 0; gout < this->ngroups(); gout++) {
      Es_ret += Es_(gin, gout);
    }
    return Es_ret;
  }

  // Operators for constructing compound cross sections
  TransportXS operator+(const TransportXS& R) const;
  TransportXS operator*(double N) const;
  TransportXS& operator+=(const TransportXS& R);
  TransportXS& operator*=(double N);

  // private:
  xt::xarray<double> Es_;    // Scattering matrix
  std::vector<double> Et_;   // Total xs
  std::vector<double> Ea_;   // Absorption xs
  std::vector<double> Ef_;   // Fission xs
  std::vector<double> nu_;   // Fission yields
  std::vector<double> chi_;  // Fission spectrum
  bool fissile_;             // Fissile indicator
};

#endif
