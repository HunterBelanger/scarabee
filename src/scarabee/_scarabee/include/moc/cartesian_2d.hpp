#ifndef CARTESIAN_2D_H
#define CARTESIAN_2D_H

#include <moc/cell.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/surface.hpp>
#include <moc/pin_cell.hpp>
#include <moc/direction.hpp>
#include <moc/vector.hpp>
#include <data/cross_section.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <variant>
#include <vector>
#include <utility>

namespace scarabee {

// Must pre-declare class so that we can have a pointer to it in the tile
class Cartesian2D;

class Cartesian2D {
 public:
  using TileFill =
      std::variant<std::shared_ptr<Cartesian2D>, std::shared_ptr<Cell>>;

  struct Tile {
    std::shared_ptr<Cartesian2D> c2d;
    std::shared_ptr<Cell> cell;

    bool valid() const {
      if (c2d && cell) return false;
      return c2d || cell;
    }

    std::size_t num_fsrs() const;

    std::size_t ngroups() const;

    std::size_t get_num_fsr_instances(std::size_t id) const;
  };

  struct TileIndex {
    std::size_t i, j;
  };

  Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
              const std::vector<std::shared_ptr<Surface>>& y_bounds);

  Cartesian2D(const std::vector<double>& dx, const std::vector<double>& dy);

  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }

  double dx() const { return x_bounds_.back()->x0() - x_bounds_.front()->x0(); }
  double dy() const { return y_bounds_.back()->y0() - y_bounds_.front()->y0(); }

  bool inside(const Vector& r, const Direction& u) const {
    return this->get_tile_index(r, u).has_value();
  }

  double distance(const Vector& r, const Direction& u) const {
    Surface x, y;
    x.type() = Surface::Type::XPlane;
    y.type() = Surface::Type::YPlane;

    x.x0() = x_min();
    const double x_min_dist = x.distance(r, u);

    x.x0() = x_max();
    const double x_max_dist = x.distance(r, u);

    y.y0() = y_min();
    const double y_min_dist = y.distance(r, u);

    y.y0() = y_max();
    const double y_max_dist = y.distance(r, u);

    return std::min(std::min(x_min_dist, x_max_dist),
                    std::min(y_min_dist, y_max_dist));
  }

  std::optional<TileIndex> get_tile_index(const Vector& r,
                                          const Direction& u) const {
    for (std::size_t i = 0; i < nx(); i++) {
      for (std::size_t j = 0; j < ny(); j++) {
        // Get the surfaces that make up our tile
        const auto& xl = x_bounds_[i];
        const auto& xh = x_bounds_[i + 1];
        const auto& yl = y_bounds_[j];
        const auto& yh = y_bounds_[j + 1];

        // Check if we are in the tile. If so, return index
        if (xl->side(r, u) == Surface::Side::Positive &&
            xh->side(r, u) == Surface::Side::Negative &&
            yl->side(r, u) == Surface::Side::Positive &&
            yh->side(r, u) == Surface::Side::Negative) {
          return TileIndex{i, j};
        }
      }
    }

    // If we get here, we weren't in any tile.
    // Return nullopt
    return std::nullopt;
  }

  Vector get_tile_center(const TileIndex& ti) const {
    // Get the lower and upper x and y bounds
    const auto& xl = x_bounds_[ti.i];
    const auto& xh = x_bounds_[ti.i + 1];
    const auto& yl = y_bounds_[ti.j];
    const auto& yh = y_bounds_[ti.j + 1];

    return Vector(0.5 * (xl->x0() + xh->x0()), 0.5 * (yl->y0() + yh->y0()));
  }

  const Tile& tile(const TileIndex& ti) const {
    check_tile_index(ti);
    return tiles_(ti.i, ti.j);
  }

  void set_tiles(const std::vector<TileFill>& fills);

  bool tiles_valid() const;

  std::shared_ptr<CrossSection> get_xs(const Vector& r,
                                       const Direction& u) const;

  UniqueFSR get_fsr(const Vector& r, const Direction& u) const;

  std::pair<UniqueFSR, Vector> get_fsr_r_local(const Vector& r,
                                               const Direction& u) const {
    auto ti = this->get_tile_index(r, u);

    if (ti.has_value() == false) {
      // We apparently are not in a tile
      return {{nullptr, 0}, r};
    }

    const auto& t = this->tile(*ti);
    const Vector tc = get_tile_center(*ti);
    Vector r_tile = r - tc;

    if (t.valid() == false) {
      std::stringstream mssg;
      mssg << "Tile for Position r = " << r << ", and Direction u = " << u
           << " is empty.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    std::pair<UniqueFSR, Vector> out{{nullptr, 0}, r_tile};

    if (t.c2d) {
      out = t.c2d->get_fsr_r_local(r_tile, u);
    } else {
      out.first = t.cell->get_fsr(r_tile, u);
    }

    // Get unique ID
    if (out.first.fsr) {
      out.first.instance +=
          fsr_offset_map_(ti->i, ti->j).find(out.first.fsr->id())->second;
    }

    return out;
  }

  std::vector<UniqueFSR> get_all_fsr_in_cell(const Vector& r,
                                             const Direction& u) const;

  std::size_t get_num_fsr_instances(std::size_t id) const;

  std::set<std::size_t> get_all_fsr_ids() const;

  void fill_fsrs(std::map<std::size_t, const FlatSourceRegion*>& fsrs) const;

  std::size_t num_fsrs() const;

  std::size_t ngroups() const { return tiles_[0].ngroups(); }

  double x_min() const { return x_bounds_.front()->x0(); }
  double x_max() const { return x_bounds_.back()->x0(); }

  double y_min() const { return y_bounds_.front()->y0(); }
  double y_max() const { return y_bounds_.back()->y0(); }

 private:
  std::vector<std::shared_ptr<Surface>> x_bounds_;
  std::vector<std::shared_ptr<Surface>> y_bounds_;
  xt::xtensor<Tile, 2> tiles_;
  xt::xtensor<std::map<std::size_t, std::size_t>, 2> fsr_offset_map_;
  std::size_t nx_, ny_;

  Cartesian2D() = default;

  void set_tile(const TileIndex& ti, const std::shared_ptr<Cartesian2D>& c2d);
  void set_tile(const TileIndex& ti, const std::shared_ptr<Cell>& cell);
  std::pair<double, double> tile_dx_dy(const TileIndex& ti) const;

  void check_tile_index(const TileIndex& ti) const {
    if (ti.i >= nx_) {
      auto mssg = "TileIndex i out of range.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }

    if (ti.j >= ny_) {
      auto mssg = "TileIndex j out of range.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Methods for getting offset map
  void make_offset_map();
};

}  // namespace scarabee

#endif
