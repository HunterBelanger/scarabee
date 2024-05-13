#ifndef SCARABEE_TRANSPORT_CROSS_SECTIONS_H
#define SCARABEE_TRANSPORT_CROSS_SECTIONS_H

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <cstdint>
#include <string>

namespace scarabee {

class TransportXS {
 public:
  TransportXS(const xt::xtensor<double, 1>& Et,
              const xt::xtensor<double, 1>& Ea,
              const xt::xtensor<double, 2>& Es,
              const xt::xtensor<double, 1>& Ef,
              const xt::xtensor<double, 1>& vEf,
              const xt::xtensor<double, 1>& chi, const std::string& name = "");

  TransportXS(const xt::xtensor<double, 1>& Et,
              const xt::xtensor<double, 1>& Ea,
              const xt::xtensor<double, 2>& Es, const std::string& name = "");

  std::size_t ngroups() const { return Et_.size(); }

  const std::string& name() const { return name_; }

  bool fissile() const { return fissile_; }

  const xt::xtensor<double, 1>& Et() const { return Et_; }

  double Et(std::size_t g) const { return Et_(g); }

  double Ea(std::size_t g) const { return Ea_(g); }

  double Ef(std::size_t g) const { return Ea_(g); }

  double vEf(std::size_t g) const { return vEf_(g); }

  double nu(std::size_t g) const { 
    const double Efg = Ef_(g);
    if (Efg > 0.) {
      return vEf_(g) / Efg;
    }
    return 0.;
  }

  double Er(std::size_t g) const { return Ea(g) + Es(g) - Es(g, g); }

  double chi(std::size_t g) const { return chi_(g); }

  double Es(std::size_t gin, std::size_t gout) const { return Es_(gin, gout); }

  double Es(std::size_t gin) const {
    return xt::sum(xt::view(Es_, gin, xt::all()))();
  }

  // Operators for constructing compound cross sections
  TransportXS operator+(const TransportXS& R) const;
  TransportXS operator*(double N) const;
  TransportXS& operator+=(const TransportXS& R);
  TransportXS& operator*=(double N);

 private:
  xt::xtensor<double, 2> Es_;   // Scattering matrix
  xt::xtensor<double, 1> Et_;   // Total xs
  xt::xtensor<double, 1> Ea_;   // Absorption xs
  xt::xtensor<double, 1> Ef_;   // Absorption xs
  xt::xtensor<double, 1> vEf_;  // Fission xs * yield
  xt::xtensor<double, 1> chi_;  // Fission spectrum
  std::string name_;
  bool fissile_;

  void check_xs();
};

}  // namespace scarabee

#endif
