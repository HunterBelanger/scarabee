#ifndef SCARABEE_NEM_DIFFUSION_DRIVER_H
#define SCARABEE_NEM_DIFFUSION_DRIVER_H

#include <diffusion_cross_section.hpp>
#include <diffusion/diffusion_geometry.hpp>

#include <Eigen/Dense>
#include <Eigen/LU>

#include <xtensor/xtensor.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <tuple>

namespace scarabee {

inline double f0(double xi) { return 1.; }

inline double f1(double xi) { return xi; }

inline double f2(double xi) { return 3. * xi * xi - 0.25; }

inline double f3(double xi) { return xi * (xi - 0.5) * (xi + 0.5); }

inline double f4(double xi) {
  return (xi * xi - 0.05) * (xi - 0.5) * (xi + 0.5);
}

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
  xt::xtensor<double, 4> flux(const xt::xtensor<double, 1>& x,
                              const xt::xtensor<double, 1>& y,
                              const xt::xtensor<double, 1>& z) const;
  xt::xtensor<double, 4> avg_flux() const;

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
  const std::size_t NG_;  // Number of groups
  const std::size_t NM_;  // Number of regions

  // Quantities required for reconstructing the flux  (kept after solution)
  xt::xtensor<double, 2> flux_avg_;  // First index is group, second is node
  xt::xtensor<double, 2> flux_x1_;
  xt::xtensor<double, 2> flux_x2_;
  xt::xtensor<double, 2> flux_y1_;
  xt::xtensor<double, 2> flux_y2_;
  xt::xtensor<double, 2> flux_z1_;
  xt::xtensor<double, 2> flux_z2_;
  xt::xtensor<Current, 2> j_outs_;
  xt::xtensor<Current, 2> j_ins_;

  // Quantites used for calculation (not needed for reconstruction)
  xt::xtensor<RMat, 2> Rmats_;  // First index is group, second is node
  xt::xtensor<PMat, 2> Pmats_;
  xt::xtensor<MomentsVector, 2> Q_;  // Source

  double keff_ = 1.;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  bool solved_{false};

  //----------------------------------------------------------------------------
  // PRIVATE METHODS
  void fill_coupling_matrices();
  void fill_source();
  void update_Jin_from_Jout(std::size_t g, std::size_t m);
  MomentsVector calc_leakage_moments(std::size_t g, std::size_t m) const;
  double calc_keff(double keff, const xt::xtensor<double, 2>& old_flux,
                   const xt::xtensor<double, 2>& new_flux) const;
  void calc_node(const std::size_t g, const std::size_t m,
                 const xt::svector<std::size_t>& geom_indx,
                 const double invs_dx, const double invs_dy,
                 const double invs_dz, const DiffusionCrossSection& xs);
  void inner_iteration();

  inline double calc_net_current(const Current& Jin, const Current& Jout,
                                 CurrentIndx indx) const {
    if (indx == CurrentIndx::XP || indx == CurrentIndx::YP ||
        indx == CurrentIndx::ZP) {
      return Jout(indx) - Jin(indx);
    } else {
      return Jin(indx) - Jout(indx);
    }
  } 

  struct P0 {
    double operator()(double /*x*/) const {
      return 1.;
    }

    double diff(double x) const {
      return 0.;
    }

    double intgr(double del) const {
      return del;
    }
  };

  struct P1 {
    double operator()(double x) const {
      return x;
    }

    double diff(double /*x*/) const {
      return 1.;
    }

    double intgr(double del) const {
      return 0.;
    }
  };

  struct P2 {
    double operator()(double x) const {
      return 0.5*(3.*x*x - 1.);
    }

    double diff(double x) const {
      return 0.5*6.*x;
    }

    double intgr(double del) const {
      return 0.;
    }
  };

  template <class Fx, class Fy>
  struct F {
    Fx fx;
    Fy fy;

    double operator()(double x, double y) const {
      return fx(x) * fy(y);
    }

    double dxiy(double dely, double x) const {
      return fx.diff(x) * fy.intgr(dely);
    }

    double dyix(double delx, double y) const {
      return fx.intgr(delx) * fy.diff(y);
    }

    double ingr(double delx, double dely) const {
      return fx.intgr(delx) * fy.intgr(dely);
    }

    double ix(double delx, double y) const {
      return fx.intgr(delx) * fy(y);
    }

    double iy(double dely, double x) const {
      return fx(x) * fy.intgr(dely);
    }
  };

  using F00 = F<P0, P0>;
  using F01 = F<P0, P1>;
  using F02 = F<P0, P2>;
  using F10 = F<P1, P0>;
  using F11 = F<P1, P1>;
  using F12 = F<P1, P2>;
  using F20 = F<P2, P0>;
  using F21 = F<P2, P1>;
  using F22 = F<P2, P2>;

  struct FluxRecon {
    std::array<double, 9> radial;
    std::array<double, 5> axial;
    double x_low, x_hi, y_low, y_hi, z_low, z_hi;

    F00 f00; F01 f01; F02 f02;
    F10 f10; F11 f11; F12 f12;
    F20 f20; F21 f21; F22 f22;

    double operator()(double x, double y, double z) const {
      const double xi_x = (x - 0.5*(x_low + x_hi)) / (0.5 * (x_hi - x_low));
      const double xi_y = (y - 0.5*(y_low + y_hi)) / (0.5 * (y_hi - y_low));
      const double xi_z = (z - 0.5*(z_low + z_hi)) / (0.5 * (z_hi - z_low));

      const double flx_z = axial[0] + axial[1]*f1(xi_z) + axial[2]*f2(xi_z) + axial[3]*f3(xi_z) + axial[4]*f4(xi_z);
      const double flx_xy = radial[0]*f00(xi_x, xi_y) + radial[1]*f01(xi_x, xi_y) + radial[2]*f02(xi_x, xi_y) + radial[3]*f10(xi_x, xi_y) + radial[4]*f11(xi_x, xi_y) + radial[5]*f12(xi_x, xi_y) + radial[6]*f20(xi_x, xi_y) + radial[7]*f21(xi_x, xi_y) + radial[8]*f22(xi_x, xi_y);
      return flx_xy * flx_z / axial[0];
    }
  };

  xt::xtensor<FluxRecon, 2> recon_params;

  FluxRecon fit_node_recon_params(std::size_t g, std::size_t m) const;

  enum class Corner {PP, PM, MP, MM};
  double eval_corner_flux(std::size_t g, std::size_t m, Corner c) const;
  double avg_corner_flux(std::size_t g, std::size_t m, Corner c) const;
};

}  // namespace scarabee

#endif
