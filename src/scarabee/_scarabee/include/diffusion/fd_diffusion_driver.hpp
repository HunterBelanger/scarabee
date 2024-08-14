#ifndef SCARABEE_FD_DIFFUSION_DRIVER_H
#define SCARABEE_FD_DIFFUSION_DRIVER_H

#include <diffusion/diffusion_data.hpp>
#include <diffusion/diffusion_geometry.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>

#include <memory>
#include <optional>
#include <tuple>

namespace scarabee {

class FDDiffusionDriver {
 public:
  FDDiffusionDriver(std::shared_ptr<DiffusionGeometry> geom);

  std::size_t ngroups() const { return geom_->ngroups(); }

  void solve();
  bool solved() const { return solved_; }

  double keff_tolerance() const { return keff_tol_; }
  void set_keff_tolerance(double ktol);

  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  double keff() const { return keff_; }

  std::tuple<xt::xarray<double>, xt::xarray<double>,
             std::optional<xt::xarray<double>>,
             std::optional<xt::xarray<double>>>
  flux() const;

  std::tuple<xt::xarray<double>, xt::xarray<double>,
             std::optional<xt::xarray<double>>,
             std::optional<xt::xarray<double>>>
  power() const;

 private:
  std::shared_ptr<DiffusionGeometry> geom_;
  xt::xtensor<double, 1> flux_;  // Flux in each MAT tile in each group
  double keff_ = 1.;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  bool solved_{false};
};

}  // namespace scarabee

#endif
