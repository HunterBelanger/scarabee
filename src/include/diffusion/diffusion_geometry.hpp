#ifndef SCARABEE_DIFFUSION_GEOMETRY_H
#define SCARABEE_DIFFUSION_GEOMETRY_H

#include <diffusion_cross_section.hpp>

#include <xtensor/xarray.hpp>

#include <cstdint>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

namespace scarabee {

class DiffusionGeometry {
 public:
  struct Tile {
    std::optional<double> albedo;
    std::shared_ptr<DiffusionCrossSection> xs;
  };

  using TileFill = std::variant<double, std::shared_ptr<DiffusionCrossSection>>;

  enum class Neighbor : std::uint8_t { XN, XP, YN, YP, ZN, ZP };

  DiffusionGeometry(const std::vector<TileFill>& tiles,
                    const std::vector<double>& dx,
                    const std::vector<std::size_t>& xdivs, double albedo_xn,
                    double albedo_xp);

  DiffusionGeometry(const std::vector<TileFill>& tiles,
                    const std::vector<double>& dx,
                    const std::vector<std::size_t>& xdivs,
                    const std::vector<double>& dy,
                    const std::vector<std::size_t>& ydivs, double albedo_xn,
                    double albedo_xp, double albedo_yn, double albedo_yp);

  std::size_t ngroups() const;
  std::size_t ndims() const {
    return tiles_.shape().size();
  }  // Number of spatial dimentions

  std::size_t nmats() const { return nmats_; }

  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }
  std::size_t nz() const { return nz_; }

  std::pair<Tile, std::optional<std::size_t>> neighbor(std::size_t m,
                                                       Neighbor n) const;
  const std::shared_ptr<DiffusionCrossSection>& mat(std::size_t m) const;
  xt::svector<std::size_t> geom_indx(std::size_t m) const;
  double volume(std::size_t m) const;

  double dx(std::size_t i) const;
  double dy(std::size_t j) const;
  double dz(std::size_t k) const;

 private:
  xt::xarray<Tile> tiles_;

  // Boundary condition tiles
  Tile xn_, xp_, yn_, yp_, zn_, zp_;

  // The bounds for each tile
  std::vector<double> tile_dx_;
  std::vector<std::size_t> x_divs_per_tile_;

  std::vector<double> tile_dy_;
  std::vector<std::size_t> y_divs_per_tile_;

  std::vector<double> tile_dz_;
  std::vector<std::size_t> z_divs_per_tile_;

  std::size_t nmats_;
  std::vector<std::size_t> mat_indx_to_flat_geom_indx_;

  std::size_t nx_, ny_, nz_;
  xt::svector<std::size_t> geom_shape_;

  xt::svector<std::size_t> geom_to_tile_indx(
      const xt::svector<std::size_t>& geo_indx) const;

  std::size_t geom_to_mat_indx(const xt::svector<std::size_t>& geo_indx) const;

  std::pair<Tile, std::optional<std::size_t>> neighbor_1d(std::size_t m,
                                                          Neighbor n) const;
  std::pair<Tile, std::optional<std::size_t>> neighbor_2d(std::size_t m,
                                                          Neighbor n) const;
  std::pair<Tile, std::optional<std::size_t>> neighbor_3d(std::size_t m,
                                                          Neighbor n) const;

  std::size_t geom_x_indx_to_tile_x_indx(std::size_t i) const;
  std::size_t geom_y_indx_to_tile_y_indx(std::size_t i) const;
  std::size_t geom_z_indx_to_tile_z_indx(std::size_t i) const;
};

}  // namespace scarabee

#endif