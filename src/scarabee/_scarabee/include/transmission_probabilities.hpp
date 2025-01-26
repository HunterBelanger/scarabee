#ifndef SCARABEE_TRANSMISSION_PROBABILITIES_H
#define SCARABEE_TRANSMISSION_PROBABILITIES_H

#include <data/cross_section.hpp>
#include <moc/boundary_condition.hpp>
#include <utils/serialization.hpp>

#include <xtensor/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>

#include <memory>

namespace scarabee {

class TransmissionProbabilities {
 public:
  TransmissionProbabilities(
      const std::vector<std::shared_ptr<CrossSection>>& xs,
      const std::vector<double>& dx, const std::vector<double>& dy,
      BoundaryCondition x_min_bc, BoundaryCondition x_max_bc,
      BoundaryCondition y_min_bc, BoundaryCondition y_max_bc);

  std::size_t ngroups() const { return ngroups_; }
  std::size_t nregions() const { return xs_.size(); }

  std::size_t nx() const { return dx_.size(); }
  std::size_t ny() const { return dy_.size(); }

  void generate_tracks(const std::uint32_t n_angles, const double d);

  void solve();
  bool solved() const { return solved_; }

  double keff_tolerance() const { return k_tol_; }
  void set_keff_tolerance(double ktol);

  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  double keff() const { return k_; }

  double flux(std::size_t i, std::size_t j, std::size_t g) const;

  xt::xtensor<double, 1> flux_spectrum(std::size_t i, std::size_t j) const;

 private:
  xt::xtensor<double, 1> dx_, dy_;
  xt::xtensor<std::shared_ptr<CrossSection>, 2> xs_;  // x, y
  xt::xtensor<double, 4> flux_incurrent_;  // Group, x, y, flux/incurrents
  xt::xtensor<double, 4> T_;               // Group, x, y, prob

  BoundaryCondition x_min_bc_, x_max_bc_, y_min_bc_, y_max_bc_;
  std::size_t ngroups_;
  double k_;
  double k_tol_;
  double flux_tol_;
  bool solved_;

  // Index definitions
  static constexpr std::size_t Flux = 0, Jxm = 1, Jxp = 2, Jym = 3, Jyp = 4;
  static constexpr std::size_t i_i = 0, xm_i = 1, xp_i = 2, ym_i = 3, yp_i = 4,
                        xm_xp = 5, xm_ym = 6, xm_yp = 7, xp_ym = 8, xp_yp = 9,
                        ym_yp = 10;

  void trace_tile(const std::size_t i, const std::size_t j,
                  const std::uint32_t NA, const double delta_phi,
                  const double d);
  void solve_tile(const std::size_t g, const std::size_t i, const std::size_t j,
                  xt::xtensor<double, 4>& flux_incurrent,
                  const xt::xtensor<double, 3>& Q) const;
  void fill_src(xt::xtensor<double, 3>& Q,
                const xt::xtensor<double, 4>& flux_incurrent);
  double calc_keff(const xt::xtensor<double, 4>& new_flux_incurrent,
                   const xt::xtensor<double, 4>& flux_incurrent) const;
  void sweep(const xt::xtensor<double, 3>& Q,
             xt::xtensor<double, 4>& flux_incurrent) const;

  friend class cereal::access;
  TransmissionProbabilities() {}
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(dx_), CEREAL_NVP(dy_), CEREAL_NVP(xs_),
        CEREAL_NVP(flux_incurrent_), CEREAL_NVP(T_), CEREAL_NVP(x_min_bc_),
        CEREAL_NVP(x_max_bc_), CEREAL_NVP(y_min_bc_), CEREAL_NVP(y_max_bc_),
        CEREAL_NVP(k_), CEREAL_NVP(k_tol_), CEREAL_NVP(flux_tol_),
        CEREAL_NVP(solved_));
  }
};

}  // namespace scarabee

#endif