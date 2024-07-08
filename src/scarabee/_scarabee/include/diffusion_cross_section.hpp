#ifndef SCARABEE_DIFFUISON_CROSS_SECTIONS_H
#define SCARABEE_DIFFUISON_CROSS_SECTIONS_H

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <cstdint>
#include <string>

namespace scarabee {

class DiffusionCrossSection {
 public:
  DiffusionCrossSection(const xt::xtensor<double, 1>& D,
                        const xt::xtensor<double, 1>& Ea,
                        const xt::xtensor<double, 2>& Es,
                        const xt::xtensor<double, 1>& Ef,
                        const xt::xtensor<double, 1>& vEf,
                        const xt::xtensor<double, 1>& chi,
                        const std::string& name = "");

  DiffusionCrossSection(const xt::xtensor<double, 1>& D,
                        const xt::xtensor<double, 1>& Ea,
                        const xt::xtensor<double, 2>& Es,
                        const std::string& name = "");

  std::size_t ngroups() const { return D_.size(); }

  const std::string& name() const { return name_; }
  void set_name(const std::string& new_name) { name_ = new_name; }

  bool fissile() const { return fissile_; }

  double D(std::size_t g) const { return D_(g); }

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

  double Er(std::size_t g) const { return Ea(g) + Es(g) - Es(g, g); }

  double chi(std::size_t g) const { return chi_(g); }

  double Es(std::size_t gin, std::size_t gout) const { return Es_(gin, gout); }

  double Es(std::size_t gin) const {
    return xt::sum(xt::view(Es_, gin, xt::all()))();
  }

 private:
  xt::xtensor<double, 2> Es_;   // Scattering matrix
  xt::xtensor<double, 1> D_;    // Diffusion coefficients
  xt::xtensor<double, 1> Ea_;   // Absorption xs
  xt::xtensor<double, 1> Ef_;   // Fission xs
  xt::xtensor<double, 1> vEf_;  // Fission xs * yield
  xt::xtensor<double, 1> chi_;  // Fission spectrum
  std::string name_;
  bool fissile_;

  void check_xs();
};

}  // namespace scarabee

#endif
