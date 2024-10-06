#ifndef MOC_DRIVER_H
#define MOC_DRIVER_H

#include <moc/cartesian_2d.hpp>
#include <moc/cmfd.hpp>
#include <moc/boundary_condition.hpp>
#include <moc/simulation_mode.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/track.hpp>
#include <moc/quadrature/polar_quadrature.hpp>
#include <utils/spherical_harmonics.hpp>

#include <xtensor/xtensor.hpp>

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
            BoundaryCondition ymax = BoundaryCondition::Reflective,
            bool anisotropic = false);

  std::shared_ptr<Cartesian2D> geometry() const { return geometry_; }

  const std::shared_ptr<CMFD>& cmfd() const { return cmfd_; }
  void set_cmfd(std::shared_ptr<CMFD> cmfd) { cmfd_ = cmfd; }

  bool drawn() const { return !angle_info_.empty(); }

  const PolarQuadrature& polar_quadrature() const { return polar_quad_; }

  std::size_t ngroups() const { return ngroups_; }

  std::size_t num_spherical_harmonics() const { return N_lj_; }

  double keff() const { return keff_; }

  SimulationMode& sim_mode() { return mode_; }
  const SimulationMode& sim_mode() const { return mode_; }

  double keff_tolerance() const { return keff_tol_; }
  void set_keff_tolerance(double ktol);

  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  void generate_tracks(std::uint32_t n_angles, double d,
                       PolarQuadrature polar_quad);

  void solve();
  bool solved() const { return solved_; }

  std::shared_ptr<CrossSection> homogenize() const;
  std::shared_ptr<CrossSection> homogenize(
      const std::vector<std::size_t>& regions) const;

  xt::xtensor<double, 1> homogenize_flux_spectrum() const;
  xt::xtensor<double, 1> homogenize_flux_spectrum(
      const std::vector<std::size_t>& regions) const;

  void apply_criticality_spectrum(const xt::xtensor<double, 1>& flux);

  std::size_t size() const;
  std::size_t nfsr() const { return this->size(); }
  std::size_t nregions() const { return this->nfsr(); }

  double flux(const Vector& r, const Direction& u, std::size_t g,
              std::size_t lj = 0) const;
  double flux(std::size_t i, std::size_t g, std::size_t lj = 0) const;

  double volume(const Vector& r, const Direction& u) const;
  double volume(std::size_t i) const;

  const std::shared_ptr<CrossSection>& xs(const Vector& r,
                                          const Direction& u) const;
  const std::shared_ptr<CrossSection>& xs(std::size_t i) const;

  UniqueFSR get_fsr(const Vector& r, const Direction& u) const;

  std::size_t get_fsr_indx(const UniqueFSR& fsr) const;

  std::vector<std::size_t> get_all_fsr_in_cell(const Vector& r,
                                               const Direction& u) const;

  std::vector<std::pair<std::size_t, double>> trace_fsr_segments(
      const Vector r_start, const Direction& u) const;

  void set_extern_src(const Vector& r, const Direction& u, std::size_t g,
                      double src);
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

  void save_hdf5(const std::string& fname, const std::string& group) const;

 private:
  struct AngleInfo {
    double phi;                  // Azimuthal angle for track
    double d;                    // Spacing for trackings of this angle
    double wgt;                  // Weight for tracks with this angle
    std::uint32_t nx;            // Number of tracks starting on the -y boundary
    std::uint32_t ny;            // Number of tracks starting on the -x boundary
    std::size_t forward_index;   // azimuthal angle index in forward
    std::size_t backward_index;  // azimuthal angle index in backward
  };

  std::vector<AngleInfo> angle_info_;       // Information for all angles
  std::vector<std::vector<Track>> tracks_;  // All tracks, indexed by angle
  std::shared_ptr<Cartesian2D> geometry_;   // Geometry for the problem
  std::shared_ptr<CMFD> cmfd_;              // CMFD for acceleration
  PolarQuadrature polar_quad_;              // Polar quadrature
  SphericalHarmonics sph_harm_;             // Spherical harmonics
  xt::xtensor<double, 3>
      flux_;  // Indexed by group, FSR, and spherical harmonic
  xt::xtensor<double, 2> extern_src_;  // Indexed by group then FSR
  std::vector<const FlatSourceRegion*> fsrs_;
  std::map<std::size_t, std::size_t> fsr_offsets_;  // Indexed by id -> offset
  std::size_t ngroups_;
  std::size_t nfsrs_;
  std::size_t n_pol_angles_;
  double flux_tol_ = 1.E-5;
  double keff_tol_ = 1.E-5;
  double keff_ = 1.;
  BoundaryCondition x_min_bc_, x_max_bc_, y_min_bc_, y_max_bc_;
  std::size_t max_L_ = 0;     // max-legendre-order in scattering moments
  std::size_t N_lj_ = 1;      // total number of j (-l ro l)
  bool anisotropic_ = false;  // to account for anisotropic scattering
  SimulationMode mode_{SimulationMode::Keff};
  bool solved_{false};

  void generate_azimuthal_quadrature(std::uint32_t n_angles, double d);
  void generate_tracks();
  void set_track_ends_bcs();
  void allocate_track_fluxes();
  void segment_renormalization();

  // isotropic
  void solve_isotropic();
  void sweep(xt::xtensor<double, 3>& flux, const xt::xtensor<double, 2>& src);
  void fill_source(xt::xtensor<double, 2>& src,
                   const xt::xtensor<double, 3>& flux) const;

  // anisotropic
  void solve_anisotropic();
  void sweep_anisotropic(xt::xtensor<double, 3>& flux,
                         const xt::xtensor<double, 3>& src);
  void fill_source_anisotropic(xt::xtensor<double, 3>& src,
                               const xt::xtensor<double, 3>& flux) const;

  double calc_keff(const xt::xtensor<double, 3>& flux,
                   const xt::xtensor<double, 3>& old_flux) const;
};

}  // namespace scarabee

#endif
