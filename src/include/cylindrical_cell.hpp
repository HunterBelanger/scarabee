#ifndef SCARABEE_CYLINDRICAL_CELL_H
#define SCARABEE_CYLINDRICAL_CELL_H

#include <cross_section.hpp>
#include <utils/constants.hpp>

#include <xtensor/xtensor.hpp>

#include <memory>
#include <vector>

namespace scarabee {

class CylindricalCell {
 public:
  CylindricalCell(const std::vector<double>& radii,
                  const std::vector<std::shared_ptr<CrossSection>>& mats);
  bool solved() const { return solved_; }
  void solve();

  std::size_t ngroups() const { return ngroups_; }
  std::size_t nregions() const { return radii_.size(); }

  // Surface area of cell
  double Sb() const { return 2. * PI * radii_.back(); }

  // Volume of region i
  double volume(std::size_t i) const { return vols_[i]; }

  // Outer radius of region i
  double radius(std::size_t i) const { return radii_[i]; }

  double Y(double a, std::uint32_t g, std::size_t i) const {
    const double Yi = Y_(g, i);
    const double Gamma = Gamma_(g);

    return Yi / (1. - a * (1. - Gamma));
  }

  double x(std::uint32_t g, std::size_t i) const {
    return 0.25 * Sb() * volume(i) * Y_(g, i);
  }

  double X(double a, std::uint32_t g, std::size_t i, std::size_t k) const {
    const double xk = x(g, k);
    const double Xik = X_(g, i, k);

    return Xik + a * xk * Y(a, g, i);
  }

  double Gamma(std::uint32_t g) const { return Gamma_(g); }

  double p(std::uint32_t g, std::size_t i, std::size_t j) const {
    return p_(g, i, j);
  }

  const std::shared_ptr<CrossSection>& xs(std::size_t i) const {
    return mats_[i];
  }

 private:
  xt::xtensor<double, 3> p_;
  xt::xtensor<double, 3> X_;
  xt::xtensor<double, 2> Y_;
  xt::xtensor<double, 1> Gamma_;  // Multicollision blackness in each group
  std::vector<double> radii_;
  std::vector<double> vols_;
  std::vector<std::shared_ptr<CrossSection>> mats_;
  std::size_t ngroups_;
  bool solved_;

  void calculate_collision_probabilities();
  double calculate_S_ij(std::size_t i, std::size_t j, std::uint32_t g) const;
  void solve_systems();
};

}  // namespace scarabee

#endif
