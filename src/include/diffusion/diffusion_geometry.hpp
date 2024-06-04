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

    enum class Neighbor : std::uint8_t {XN, XP, YN, YP, ZN, ZP};

    DiffusionGeometry(const std::vector<double>& dx, double albedo_xn, double albedo_xp);

    void set_tiles(const std::vector<TileFill>& tiles);

    std::size_t ngroups() const;
    std::size_t ndims() const { return tiles_.shape().size(); } // Number of spatial dimentions

    std::size_t ntiles() const { return tiles_.size(); }
    std::size_t nmats() const { return n_mat_tiles_; }

    std::size_t nx() const { return x_bounds_.size() > 0 ? x_bounds_.size() - 1 : 0; }
    std::size_t ny() const { return y_bounds_.size() > 0 ? y_bounds_.size() - 1 : 0; }
    std::size_t nz() const { return z_bounds_.size() > 0 ? z_bounds_.size() - 1 : 0; }

    std::pair<Tile, std::optional<std::size_t>> neighbor(std::size_t m, Neighbor n) const;
    const std::shared_ptr<DiffusionCrossSection>& mat(std::size_t m) const;
    std::vector<std::size_t> mat_indxs(std::size_t m) const;
    double volume(std::size_t m) const;

    double dx(std::size_t i) const { return x_bounds_[i+1] - x_bounds_[i]; }
    double dy(std::size_t j) const { return y_bounds_[j+1] - y_bounds_[j]; }
    double dz(std::size_t k) const { return z_bounds_[k+1] - z_bounds_[k]; }

  private:
    xt::xarray<Tile> tiles_;

    // Boundary condition tiles
    Tile xn_, xp_, yn_, yp_, zn_, zp_;

    // The bounds for each tile
    std::vector<double> x_bounds_;
    std::vector<double> y_bounds_;
    std::vector<double> z_bounds_;

    std::size_t n_mat_tiles_;

    std::vector<std::size_t> mat_indx_to_flat_tile_indx_;

    void set_tiles_1d(const std::vector<TileFill>& tiles);
};

}

#endif