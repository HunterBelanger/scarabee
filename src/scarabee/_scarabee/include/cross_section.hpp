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
               const xt::xtensor<double, 1>& Dtr,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 3>& Es,
               const xt::xtensor<double, 1>& Ef,
               const xt::xtensor<double, 1>& vEf,
               const xt::xtensor<double, 1>& chi, const std::string& name = "");

  CrossSection(const xt::xtensor<double, 1>& Etr,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 2>& Es_tr,
               const std::string& name = "");

  CrossSection(const xt::xtensor<double, 1>& Et,
               const xt::xtensor<double, 1>& Dtr,
               const xt::xtensor<double, 1>& Ea,
               const xt::xtensor<double, 3>& Es, const std::string& name = "");

  std::size_t ngroups() const { return Etr_.size(); }

  const std::string& name() const { return name_; }
  void set_name(const std::string& new_name) { name_ = new_name; }

  bool fissile() const { return fissile_; }

  bool anisotropic() const { return Es_.shape()[0] > 1; }

  std::size_t max_legendre_order() const { return Es_.shape()[0]-1; }

  const xt::xtensor<double, 1>& Etr() const { return Etr_; }

  double Etr(std::size_t g) const { return Etr_(g); }

  double Dtr(std::size_t g) const { return Dtr_(g); }

  double Et(std::size_t g) const {
    return Etr_(g) + Dtr_(g);
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
    return Es_(0, gin, gout);
  }

  double Es_tr(std::size_t gin) const {
    return xt::sum(xt::view(Es_, 0, gin, xt::all()))();
  }

  double Es(std::size_t l, std::size_t gin, std::size_t gout) const {
    if (l > this->max_legendre_order())
      return 0.;

    double Es_l_gin_gout = Es_(l, gin, gout);
    if (l == 0 && gin == gout) {
      Es_l_gin_gout += Dtr_(gin);
    }
    return Es_l_gin_gout;
  }

  double Es(std::size_t l, std::size_t gin) const {
    if (l > this->max_legendre_order())
      return 0.;

    double Es = xt::sum(xt::view(Es_, l, gin, xt::all()))();
    if (l == 0) {
      Es += Dtr_[gin];
    }
    return Es;
  }

  std::shared_ptr<CrossSection> condense(
      const std::vector<std::pair<std::size_t, std::size_t>>& groups,
      const xt::xtensor<double, 1>& flux) const;

  // Operators for constructing compound cross sections
  CrossSection operator+(const CrossSection& R) const;
  CrossSection operator*(double N) const;
  CrossSection& operator+=(const CrossSection& R);
  CrossSection& operator*=(double N);

 private:
  xt::xtensor<double, 1> Etr_; // Transport xs
  xt::xtensor<double, 1> Dtr_; // Transport Correction xs
  xt::xtensor<double, 1> Ea_;  // Absorption xs
  xt::xtensor<double, 1> Ef_;  // Fission xs
  xt::xtensor<double, 1> vEf_; // Fission xs * yield
  xt::xtensor<double, 1> chi_; // Fission spectrum
  // Scattering Matrices. [0,:,:] is the TRANSPORT corrected P0 scatter matrix
  // while [1,:,:] is the P1 scatter matrix, [2,:,:] is the P2 matrix, etc.
  xt::xtensor<double, 3> Es_;
  std::string name_;
  bool fissile_;

  void check_xs();
};

}  // namespace scarabee

#endif
