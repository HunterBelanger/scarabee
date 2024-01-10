#ifndef SCARABEE_CYLINDRICAL_CELL_H
#define SCARABEE_CYLINDRICAL_CELL_H

#include <mg_cross_sections.hpp>
#include <utils/constants.hpp>

#include <ndarray.hpp>

#include <vector>

class CylindricalCell {
 public:
  CylindricalCell(const std::vector<double>& radii,
                  const std::vector<std::shared_ptr<MGCrossSections>>& mats);
  bool solved() const { return solved_; }
  void solve();

  std::uint32_t ngroups() const { return ngroups_; }
  std::size_t nregions() const { return radii_.size(); }

  // Surface area of cell
  double S() const { return PI * radii_.back() * radii_.back(); }

  // Volume of region i
  double V(std::size_t i) const { return vols_[i]; }

  // Outer radius of region i
  double R(std::size_t i) const { return radii_[i]; }

  double Y(double a, std::uint32_t g, std::size_t i) const {
    const double Yi = Y_(g, i);
    const double Gamma = Gamma_[g];

    return Yi / (1. - a * (1 - Gamma));
  }

  double X(double a, std::uint32_t g, std::size_t i, std::size_t k) const {
    const double xk = 0.25 * S() * vols_[k] * Y_(g, k);
    const double Xik = X_(g, i, k);
    const double Yi = Y_(g, i);
    const double Gamma = Gamma_[g];

    return Xik + (a * xk * (Yi / (1. - a * (1. - Gamma))));
  }

  const MGCrossSections& mat(std::size_t i) const { return *mats_[i]; }

  void print_p() const;

 private:
  NDArray<double> p_;
  NDArray<double> X_;
  NDArray<double> Y_;
  std::vector<double> Gamma_;  // Multicollision blackness in each group
  std::vector<double> radii_;
  std::vector<double> vols_;
  std::vector<std::shared_ptr<MGCrossSections>> mats_;
  std::uint32_t ngroups_;
  bool solved_;

  void calculate_collision_probabilities();
  double calculate_S_ij(std::size_t i, std::size_t j, std::uint32_t g) const;
  void solve_systems();
};

#endif
