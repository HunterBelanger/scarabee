#ifndef SCARABEE_FD_DIFFUSION_DRIVER_H
#define SCARABEE_FD_DIFFUSION_DRIVER_H

#include <diffusion/diffusion_data.hpp>
#include <diffusion/diffusion_geometry.hpp>
#include <utils/serialization.hpp>

#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>

#include <memory>
#include <optional>
#include <tuple>

namespace scarabee {

class FDDiffusionDriver {
 public:
  FDDiffusionDriver(std::shared_ptr<DiffusionGeometry> geom);

  std::shared_ptr<DiffusionGeometry> geometry() const { return geom_; }

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

  void save(const std::string& fname);
  static std::unique_ptr<FDDiffusionDriver> load(const std::string& fname);

 private:
  std::shared_ptr<DiffusionGeometry> geom_;
  xt::xtensor<double, 1> flux_;  // Flux in each MAT tile in each group
  double keff_ = 1.;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  bool solved_{false};

  friend class cereal::access;
  FDDiffusionDriver() {}
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(geom_), CEREAL_NVP(flux_), CEREAL_NVP(keff_),
        CEREAL_NVP(flux_tol_), CEREAL_NVP(keff_tol_), CEREAL_NVP(solved_));
  }
};

}  // namespace scarabee

#endif
