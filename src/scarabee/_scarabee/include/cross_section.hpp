#ifndef SCARABEE_CROSS_SECTIONS_H
#define SCARABEE_CROSS_SECTIONS_H

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <cstdint>
#include <string>
#include <memory>
#include <utility>

namespace scarabee {

class CrossSection {
 public:
  CrossSection(const xt::xtensor<double, 1>& Etr,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 2>& Es_tr,
               const xt::xtensor<double, 1>& Ef,
               const xt::xtensor<double, 1>& vEf,
               const xt::xtensor<double, 1>& chi, const std::string& name = "");

  CrossSection(const xt::xtensor<double, 1>& Et,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 2>& Es,
               const xt::xtensor<double, 2>& Es1,
               const xt::xtensor<double, 1>& Ef,
               const xt::xtensor<double, 1>& vEf,
               const xt::xtensor<double, 1>& chi, const std::string& name = "");

  CrossSection(const xt::xtensor<double, 1>& Etr,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 2>& Es_tr,
               const std::string& name = "");

  CrossSection(const xt::xtensor<double, 1>& Et,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 2>& Es,
               const xt::xtensor<double, 2>& Es1, const std::string& name = "");

  std::size_t ngroups() const { return Etr_.size(); }

  const std::string& name() const { return name_; }
  void set_name(const std::string& new_name) { name_ = new_name; }

  bool fissile() const { return fissile_; }

  bool anisotropic() const { return Es1_.size() > 0; }

  const xt::xtensor<double, 1>& Etr() const { return Etr_; }

  double Etr(std::size_t g) const { return Etr_(g); }

  double Et(std::size_t g) const {
    double Et = Etr_(g);
    if (anisotropic()) {
      Et += Es1_(g, g);
    }
    return Et;
  }

  double Ea(std::size_t g) const { return Ea_(g); }

  double Ef(std::size_t g) const { return Ef_(g); }

  double vEf(std::size_t g) const { return vEf_(g); }

  double nu(std::size_t g) const {
    const double Efg = Ef_(g);
    if (Efg > 0.) {
      return vEf_(g) / Efg;
    }
    return 0.;
  }

  double Er(std::size_t g) const { return Ea(g) + Es_tr(g) - Es_tr(g, g); }

  double chi(std::size_t g) const { return chi_(g); }

  double Es_tr(std::size_t gin, std::size_t gout) const {
    return Es_tr_(gin, gout);
  }

  double Es1(std::size_t gin, std::size_t gout) const {
    if (anisotropic()) return Es1_(gin, gout);
    return 0.;
  }

  double Es(std::size_t gin, std::size_t gout) const {
    double Es_gin_gout = Es_tr_(gin, gout);
    if (anisotropic() && gin == gout) {
      Es_gin_gout += Es1_(gin, gout);
    }
    return Es_gin_gout;
  }

  double Es_tr(std::size_t gin) const {
    return xt::sum(xt::view(Es_tr_, gin, xt::all()))();
  }

  double Es(std::size_t gin) const {
    double Es = xt::sum(xt::view(Es_tr_, gin, xt::all()))();
    if (anisotropic()) {
      Es += Es1_(gin, gin);
    }
    return Es;
  }

  std::shared_ptr<CrossSection> condense(
      const std::vector<std::pair<std::size_t, std::size_t>>& groups,
      const std::vector<double>& flux) const;

  // Operators for constructing compound cross sections
  CrossSection operator+(const CrossSection& R) const;
  CrossSection operator*(double N) const;
  CrossSection& operator+=(const CrossSection& R);
  CrossSection& operator*=(double N);

 private:
  xt::xtensor<double, 2> Es_tr_;  // Transport Corrected Scattering matrix
  xt::xtensor<double, 2> Es1_;    // P1 Scattering matrix
  xt::xtensor<double, 1> Etr_;    // Transport xs
  xt::xtensor<double, 1> Ea_;     // Absorption xs
  xt::xtensor<double, 1> Ef_;     // Fission xs
  xt::xtensor<double, 1> vEf_;    // Fission xs * yield
  xt::xtensor<double, 1> chi_;    // Fission spectrum
  std::string name_;
  bool fissile_;

  void check_xs();
};

}  // namespace scarabee

#endif
