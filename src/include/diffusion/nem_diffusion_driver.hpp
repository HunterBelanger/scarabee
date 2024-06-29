#ifndef SCARABEE_NEM_DIFFUSION_DRIVER_H
#define SCARABEE_NEM_DIFFUSION_DRIVER_H

#include <diffusion_cross_section.hpp>
#include <diffusion/diffusion_geometry.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>

#include <memory>
#include <tuple>

namespace scarabee {

class NEMDiffusionDriver {
 public:
  NEMDiffusionDriver(std::shared_ptr<DiffusionGeometry> geom);

  std::size_t ngroups() const { return geom_->ngroups(); }

  void solve();
  bool solved() const { return solved_; }

  double keff_tolerance() const { return keff_tol_; }
  void set_keff_tolerance(double ktol);

  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  double keff() const { return keff_; }

  double flux(double x, double y, double z, std::size_t g) const;
  xt::xarray<double> flux(xt::xarray<double> x, xt::xarray<double> y, xt::xarray<double> z) const;
  
  double power(double x, double y, double z) const;
  xt::xarray<double> power(xt::xarray<double> x, xt::xarray<double> y, xt::xarray<double> z) const;

 private:
  //----------------------------------------------------------------------------
  // TYPES
  // Container for node averaged flux and moments
  struct Flux {
    double avg, x1, x2, y1, y2, z1, z2;
  };
  // Definition of coupling matrix types
  using Current = xt::xtensor_fixed<double, xt::xshape<6, 1>>;
  enum CurrentIndx : ::std::size_t {XP = 0, XM = 1, YP = 2, YM = 3, ZP = 4, ZM = 5};

  using MomentsVector = xt::xtensor_fixed<double, xt::xshape<7, 1>>;
  enum MomentIndx : ::std::size_t {AVG = 0, X1 = 1, Y1 = 2, Z1 = 3, X2 = 4, Y2 = 5, Z2 = 6};

  using RMat = xt::xtensor_fixed<double, xt::xshape<6,6>>;
  using PMat = xt::xtensor_fixed<double, xt::xshape<6,7>>;
  
  //----------------------------------------------------------------------------
  // PRIVATE MEMBERS
  std::shared_ptr<DiffusionGeometry> geom_;

  // Quantities required for reconstructing the flux  (kept after solution)
  xt::xtensor<Flux, 2> flux_;      // First index is group, second is node
  xt::xtensor<Current, 2> j_outs_; // First index is group, second is node
  xt::xtensor<Current, 2> j_ins_;  // First index is group, second is node

  // Quantites used for calculation (not needed for reconstruction)
  xt::xtensor<RMat, 2> Rmats_; // First index is group, second is node
  xt::xtensor<PMat, 2> Pmats_; // First index is group, second is node

  double keff_ = 1.;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  bool solved_{false};
};

}  // namespace scarabee

#endif
