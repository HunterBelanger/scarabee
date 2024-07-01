#ifndef SCARABEE_NEM_DIFFUSION_DRIVER_H
#define SCARABEE_NEM_DIFFUSION_DRIVER_H

#include <diffusion_cross_section.hpp>
#include <diffusion/diffusion_geometry.hpp>

#include <Eigen/Dense>

#include <xtensor/xtensor.hpp>

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
  xt::xarray<double> flux(xt::xarray<double> x, xt::xarray<double> y,
                          xt::xarray<double> z) const;

  double power(double x, double y, double z) const;
  xt::xarray<double> power(xt::xarray<double> x, xt::xarray<double> y,
                           xt::xarray<double> z) const;

 private:
  //----------------------------------------------------------------------------
  // TYPES
  // Definition of coupling matrix types
  using Current = Eigen::Matrix<double, 6, 1>;
  enum CurrentIndx : ::std::size_t {
    XP = 0,
    XM = 1,
    YP = 2,
    YM = 3,
    ZP = 4,
    ZM = 5
  };

  using MomentsVector = Eigen::Matrix<double, 7, 1>;
  enum MomentIndx : ::std::size_t {
    AVG = 0,
    X1 = 1,
    Y1 = 2,
    Z1 = 3,
    X2 = 4,
    Y2 = 5,
    Z2 = 6
  };

  using RMat = Eigen::Matrix<double, 6, 6>;
  using PMat = Eigen::Matrix<double, 6, 7>;

  //----------------------------------------------------------------------------
  // PRIVATE MEMBERS
  std::shared_ptr<DiffusionGeometry> geom_;
  const std::size_t NG_; // Number of groups
  const std::size_t NM_; // Number of regions

  // Quantities required for reconstructing the flux  (kept after solution)
  xt::xtensor<double,  2> flux_avg_; // First index is group, second is node
  xt::xtensor<double,  2> flux_x1_;
  xt::xtensor<double,  2> flux_x2_;
  xt::xtensor<double,  2> flux_y1_;
  xt::xtensor<double,  2> flux_y2_;
  xt::xtensor<double,  2> flux_z1_;
  xt::xtensor<double,  2> flux_z2_;
  xt::xtensor<Current, 2> j_outs_;
  xt::xtensor<Current, 2> j_ins_;

  // Quantites used for calculation (not needed for reconstruction)
  xt::xtensor<RMat, 2> Rmats_;  // First index is group, second is node
  xt::xtensor<PMat, 2> Pmats_;
  xt::xtensor<MomentsVector, 2> Q_; // Source

  double keff_ = 1.;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  bool solved_{false};

  //----------------------------------------------------------------------------
  // PRIVATE METHODS
  void fill_coupling_matrices();
  void fill_source();
  static double calc_net_current(const Current& Jin, const Current& Jout, CurrentIndx indx);
  void update_Jin_from_Jout(std::size_t g, std::size_t m);
  MomentsVector calc_leakage_moments(std::size_t g, std::size_t m) const;
  double calc_keff(double keff, const xt::xtensor<double, 2>& old_flux, const xt::xtensor<double, 2>& new_flux) const;
  double calc_node(const std::size_t g, const std::size_t m, const xt::svector<std::size_t>& geom_indx, const double invs_dx, const double invs_dy, const double invs_dz, const DiffusionCrossSection& xs);
  void inner_iteration(xt::xtensor<double, 2>& errors);
};

}  // namespace scarabee

#endif
