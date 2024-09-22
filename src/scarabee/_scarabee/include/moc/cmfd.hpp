#ifndef CMFD_H
#define CMFD_H

#include <moc/surface.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>

#include <xtensor/xtensor.hpp>

#include <array>
#include <utility>
#include <optional>
#include <vector>
#include <set>

namespace scarabee {

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
  std::optional<std::size_t> get_surface(const Vector& r,
                                         const Direction& u) const;

  void insert_fsr(const std::array<std::size_t, 2>& tile, std::size_t fsr);
  void pack_fsr_lists();

  std::size_t moc_to_cmfd_group(std::size_t g) const;

  double& current(const std::size_t G, const std::size_t surface);
  const double& current(const std::size_t G, const std::size_t surface) const;

  void tally_current(double aflx, const Direction& u, std::size_t G, const std::size_t surf);

  void zero_currents() { surface_currents_.fill(0.); }

 private:
  std::vector<Surface> x_bounds_;
  std::vector<Surface> y_bounds_;
  std::vector<std::size_t> moc_to_cmfd_group_map_;
  std::size_t nx_, ny_;
  std::size_t nx_surfs_, ny_surfs_;

  // List of flat source region indices for each CMFD cell
  std::vector<std::set<std::size_t>> temp_fsrs_;
  std::vector<std::vector<std::size_t>> fsrs_;

  // This contains the net current for every possible surface in every group.
  // The first index is the CMFG group, and the second is the surface ID.
  // Surfaces are ordered as all x surfaces, then all y surfaces.
  // Number of surfaces is then ny_*x_bounds_.size() + nx_*y_bounds_.size().
  xt::xtensor<double, 2> surface_currents_;  // group, surface
};

}  // namespace scarabee

#endif