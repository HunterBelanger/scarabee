#ifndef SCARABEE_CMFD_H
#define SCARABEE_CMFD_H

#include <moc/surface.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>
#include <data/diffusion_cross_section.hpp>

#include <xtensor/xtensor.hpp>
#include <Eigen/Sparse>

#include <array>
#include <memory>
#include <utility>
#include <optional>
#include <vector>
#include <set>

namespace scarabee {

class MOCDriver;

struct CMFDSurfaceCrossing {
//left, right, bottom, top, top right, bottom right, bottom left, top left 

enum class Type : std::uint8_t {XN, XP, YN, YP, TP, BR, BL, TL};
std::size_t cell_index_;
bool is_valid_;
Type crossing_;

constexpr explicit operator bool() const noexcept { return is_valid; }

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

  std::optional<std::size_t> get_surface(const Vector& r,
                                         const Direction& u) const;

  void insert_fsr(const std::array<std::size_t, 2>& tile, std::size_t fsr);
  void insert_fsr(std::size_t tile_indx, std::size_t fsr);
  void pack_fsr_lists();

  std::size_t moc_to_cmfd_group(std::size_t g) const;

  double& current(const std::size_t G, const std::size_t surface);
  const double& current(const std::size_t G, const std::size_t surface) const;

  void tally_current(double aflx, const Direction& u, std::size_t G,
                     const std::size_t surf);

  void zero_currents() {
    surface_currents_.fill(0.);
    surface_currents_normalized_ = false;
  }

  void solve(MOCDriver& moc, double keff);

 private:
  std::vector<double> dx_, dy_;
  std::vector<Surface> x_bounds_;
  std::vector<Surface> y_bounds_;
  std::vector<std::size_t> moc_to_cmfd_group_map_;
  std::vector<std::pair<std::size_t, std::size_t>> group_condensation_;
  std::size_t nx_, ny_, ng_;
  std::size_t nx_surfs_, ny_surfs_;

  // List of flat source region indices for each CMFD cell
  std::vector<std::set<std::size_t>> temp_fsrs_;
  std::vector<std::vector<std::size_t>> fsrs_;

  // This contains the net current for every possible surface in every group.
  // The first index is the CMFG group, and the second is the surface ID.
  // Surfaces are ordered as all x surfaces, then all y surfaces.
  // Number of surfaces is then ny_*x_bounds_.size() + nx_*y_bounds_.size().
  xt::xtensor<double, 2> surface_currents_;  // group, surface
  bool surface_currents_normalized_ = false;

  xt::xtensor<std::shared_ptr<DiffusionCrossSection>, 2> xs_;
  xt::xtensor<double, 3> Et_; // g, i, j
  xt::xtensor<double, 3> D_transp_corr_; // g, i, j
  xt::xtensor<double, 3> flux_;  // g, x, y

  Eigen::SparseMatrix<double> M;  // Loss Matrix

  void normalize_currents();
  void compute_homogenized_xs_and_flux(const MOCDriver& moc);
  void check_neutron_balance(const std::size_t i, const std::size_t j, std::size_t g, const double keff) const;
};

}  // namespace scarabee

#endif