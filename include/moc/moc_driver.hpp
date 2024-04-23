#ifndef MOC_DRIVER_H
#define MOC_DRIVER_H

#include <moc/cartesian_2d.hpp>
#include <moc/boundary_condition.hpp>
#include <moc/simulation_mode.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/track.hpp>
#include <moc/quadrature/polar_quadrature.hpp>

#include <map>
#include <memory>
#include <vector>

namespace scarabee {

class MOCDriver {
 public:
  MOCDriver(std::shared_ptr<Cartesian2D> geometry,
            BoundaryCondition xmin = BoundaryCondition::Reflective,
            BoundaryCondition xmax = BoundaryCondition::Reflective,
            BoundaryCondition ymin = BoundaryCondition::Reflective,
            BoundaryCondition ymax = BoundaryCondition::Reflective);

  bool drawn() const { return !angle_info_.empty(); }

  const PolarQuadrature& polar_quadrature() const { return polar_quad_; }

  std::size_t ngroups() const { return ngroups_; }

  double keff() const { return keff_; }

  SimulationMode& sim_mode() { return mode_; }
  const SimulationMode& sim_mode() const { return mode_; }

  double keff_tolerance() const { return keff_tol_; }
  void set_keff_tolerance(double ktol);

  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  void generate_tracks(std::uint32_t n_angles, double d,
                       PolarQuadrature polar_quad, bool precalc_exps = true);

  void solve();
  bool solved() const { return solved_; }

  double get_flux(std::size_t g, const Vector& r, const Direction& u) const;

  std::shared_ptr<TransportXS> get_xs(const Vector& r,
                                      const Direction& u) const;

  UniqueFSR get_fsr(const Vector& r, const Direction& u) const;

  std::size_t get_fsr_indx(const UniqueFSR& fsr) const;

  void set_extern_src(const Vector& r, const Direction& u, std::size_t g, double src);
  double extern_src(const Vector& r, const Direction& u, std::size_t g) const;

  void set_extern_src(std::size_t i, std::size_t g, double src);
  double extern_src(std::size_t i, std::size_t g) const;

  BoundaryCondition& x_min_bc() { return x_min_bc_; }
  const BoundaryCondition& x_min_bc() const { return x_min_bc_; }

  BoundaryCondition& x_max_bc() { return x_max_bc_; }
  const BoundaryCondition& x_max_bc() const { return x_max_bc_; }

  BoundaryCondition& y_min_bc() { return y_min_bc_; }
  const BoundaryCondition& y_min_bc() const { return y_min_bc_; }

  BoundaryCondition& y_max_bc() { return y_max_bc_; }
  const BoundaryCondition& y_max_bc() const { return y_max_bc_; }

  double x_min() const { return geometry_->x_min(); }
  double x_max() const { return geometry_->x_max(); }
  double y_min() const { return geometry_->y_min(); }
  double y_max() const { return geometry_->y_max(); }

 private:
  struct AngleInfo {
    double phi;        // Azimuthal angle for track
    double d;          // Spacing for trackings of this angle
    double wgt;        // Weight for tracks with this angle
    std::uint32_t nx;  // Number of tracks starting on the -y boundary
    std::uint32_t ny;  // Number of tracks starting on the -x boundary
  };

  std::vector<AngleInfo> angle_info_;       // Information for all angles
  std::vector<std::vector<Track>> tracks_;  // All tracks, indexed by angle
  std::shared_ptr<Cartesian2D> geometry_;   // Geometry for the problem
  PolarQuadrature polar_quad_;              // Polar quadrature
  xt::xtensor<double, 2> flux_;             // Indexed by group then FSR
  xt::xtensor<double, 2> extern_src_;       // Indexed by group then FSR
  std::map<std::size_t, std::size_t> fsr_offsets_;
  std::size_t ngroups_;
  std::size_t nfsrs_;
  std::size_t n_pol_angles_;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  double keff_ = 1.;
  BoundaryCondition x_min_bc_, x_max_bc_, y_min_bc_, y_max_bc_;
  SimulationMode mode_ {SimulationMode::Keff};
  bool precalculated_exps_{false};
  bool solved_{false};

  void generate_azimuthal_quadrature(std::uint32_t n_angles, double d);
  void generate_tracks();
  void set_track_ends_bcs();
  void allocate_track_fluxes();
  void segment_renormalization();
  void calculate_segment_exps();

  void sweep(xt::xtensor<double, 2>& flux, const xt::xtensor<double, 2>& src);

  double calc_keff(const xt::xtensor<double, 2>& flux,
                   const xt::xtensor<double, 2>& old_flux) const;
  void fill_source(xt::xtensor<double, 2>& src,
                   const xt::xtensor<double, 2>& flux) const;
};

}  // namespace scarabee

#endif
