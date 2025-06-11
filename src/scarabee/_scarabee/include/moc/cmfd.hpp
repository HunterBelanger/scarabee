#ifndef SCARABEE_CMFD_H
#define SCARABEE_CMFD_H

#include <moc/surface.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>
#include <moc/boundary_condition.hpp>
#include <data/diffusion_cross_section.hpp>
#include <utils/simulation_mode.hpp>
#include <utils/openmp_mutex.hpp>


#include <xtensor/xtensor.hpp>
#include <Eigen/Sparse>

#include <array>
#include <memory>
#include <utility>
#include <optional>
#include <variant>
#include <vector>
#include <set>

namespace scarabee {

class MOCDriver;

struct CMFDSurfaceCrossing {
  // left, right, bottom, top, top right, bottom right, bottom left, top left

  enum class Type : std::uint8_t { XN, XP, YN, YP, TR, BR, BL, TL };
  std::size_t cell_index;
  bool is_valid;
  Type crossing;

  constexpr explicit operator bool() const noexcept { return is_valid; }

  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(cell_index), CEREAL_NVP(is_valid), CEREAL_NVP(crossing));
  }
};

class CMFD {
 public:
  CMFD(const std::vector<double>& dx, const std::vector<double>& dy,
       const std::vector<std::pair<std::size_t, std::size_t>>& groups);

  std::size_t nx() const { return nx_; };
  std::size_t ny() const { return ny_; };

  std::size_t nx_surfs() const { return nx_surfs_; }
  std::size_t ny_surfs() const { return ny_surfs_; }

  std::optional<std::array<std::size_t, 2>> get_tile(const Vector& r,
                                                     const Direction& u) const;
  std::size_t tile_to_indx(const std::array<std::size_t, 2>& tile) const;
  std::size_t tile_to_indx(const std::size_t& i, const std::size_t& j) const;

  std::array<std::size_t, 2> indx_to_tile(std::size_t cell_index);

  CMFDSurfaceCrossing get_surface(const Vector& r, const Direction& u) const;

  std::size_t get_x_neg_surf(const std::size_t i, const std::size_t j) const {
    return (nx_ + 1) * j + i;
  }
  std::size_t get_x_pos_surf(const std::size_t i, const std::size_t j) const {
    return get_x_neg_surf(i, j) + 1;
  }
  std::size_t get_y_neg_surf(const std::size_t i, const std::size_t j) const {
    return nx_surfs_ + (ny_ + 1) * i + j;
  }
  std::size_t get_y_pos_surf(const std::size_t i, const std::size_t j) const {
    return get_y_neg_surf(i, j) + 1;
  }

  void insert_fsr(const std::array<std::size_t, 2>& tile, std::size_t fsr);
  void insert_fsr(std::size_t tile_indx, std::size_t fsr);
  void pack_fsr_lists();

  std::size_t moc_to_cmfd_group(std::size_t g) const;

  double& current(const std::size_t G, const std::size_t surface);
  const double& current(const std::size_t G, const std::size_t surface) const;

  void tally_current(double aflx, const Direction& u, std::size_t G,
                     const CMFDSurfaceCrossing& surf);

  void zero_currents() {
    surface_currents_.fill(0.);
    surface_currents_normalized_ = false;
  }

  enum class TileSurf : std::uint8_t { XN, XP, YN, YP };

  void solve(MOCDriver& moc, double keff, std::size_t moc_iteration);

  double keff_tolerance() const { return keff_tol_; }
  void set_keff_tolerance(double ktol);
  double flux_tolerance() const { return flux_tol_; }
  void set_flux_tolerance(double ftol);

  void set_damping(double wd);
  void set_flux_limiting(bool user_pref) { flux_limiting_ = user_pref; }
  void set_larsen_correction(bool user_pref) { larsen_correction_ = user_pref; }
  void set_skip_moc_iterations(int num_iter) { skip_moc_iterations_ = num_iter; }

  void homogenize_ext_src(const MOCDriver& moc);

  const double& flux(const std::size_t i, const std::size_t j, const std::size_t g) const;
  double keff() const { return keff_; }

  const double& solve_time() const { return solve_time_; }
  bool solved() const { return solved_; }

 private:
  std::vector<double> dx_, dy_;
  std::vector<Surface> x_bounds_;
  std::vector<Surface> y_bounds_;
  std::vector<std::size_t> moc_to_cmfd_group_map_;
  std::vector<std::pair<std::size_t, std::size_t>> group_condensation_;
  std::size_t nx_, ny_, ng_;
  std::size_t nx_surfs_, ny_surfs_;

  bool flux_limiting_ = true;
  bool larsen_correction_ = false;
  double keff_tol_ = 1E-5;
  double flux_tol_ = 1E-5;
  double damping_ = 0.7;
  std::size_t skip_moc_iterations_ = 0;
  std::size_t moc_iteration_;
  double keff_ = 1.0;
  double solve_time_ = 0.0;
  bool solved_ = false;
  SimulationMode mode_{SimulationMode::Keff};

  // List of flat source region indices for each CMFD cell
  std::vector<std::set<std::size_t>> temp_fsrs_;
  std::vector<std::vector<std::size_t>> fsrs_;

  // This contains the net current for every possible surface in every group.
  // The first index is the CMFG group, and the second is the surface ID.
  // Surfaces are ordered as all x surfaces, then all y surfaces.
  // Number of surfaces is then ny_*x_bounds_.size() + nx_*y_bounds_.size().
  xt::xtensor<double, 2> surface_currents_;  // group, surface
  bool surface_currents_normalized_ = false;

  xt::xtensor<OpenMPMutex, 2> surface_currents_locks_;

  xt::xtensor<std::shared_ptr<DiffusionCrossSection>, 2> xs_;
  xt::xtensor<double, 3> Et_;             // g, i, j
  xt::xtensor<double, 3> flux_;           // g, x, y
  xt::xtensor<double, 2> D_transp_corr_;  // g, surf

  Eigen::VectorXd flux_cmfd_;       // g*nx_*ny_
  Eigen::VectorXd update_ratios_;   // g*nx_*ny_
  Eigen::VectorXd volumes_;         // nx_*ny_

  Eigen::SparseMatrix<double> M_;   // Loss Matrix
  Eigen::SparseMatrix<double> QM_;  // Source Matrix

  Eigen::VectorXd extern_src_; // g*nx_*ny_

  void larsen_correction(double& D, const double dx, const MOCDriver& moc) const;
  std::pair<double, double> calc_surf_diffusion_coeffs(
      std::size_t i, std::size_t j, std::size_t g, TileSurf surf,
      const MOCDriver& moc) const;
  std::variant<std::array<std::size_t, 2>, BoundaryCondition>
  find_next_cell_or_bc(std::size_t i, std::size_t j, TileSurf surf,
                       const MOCDriver& moc) const;
  double get_cmfd_tile_width(std::size_t i, std::size_t j, TileSurf surf) const;
  double get_current(std::size_t i, std::size_t j, std::size_t g,
                     TileSurf surf) const;
  void create_loss_matrix(const MOCDriver& moc);
  void create_source_matrix();
  void power_iteration(double keff);
  void fixed_source_solve();
  void update_moc_fluxes(MOCDriver& moc);
  void normalize_currents();
  void compute_homogenized_xs_and_flux(const MOCDriver& moc);
  void check_neutron_balance(const std::size_t i, const std::size_t j,
                             std::size_t g, const double keff) const;
};

}  // namespace scarabee

#endif