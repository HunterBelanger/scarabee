#ifndef CARTESIAN_2D_H
#define CARTESIAN_2D_H

#include <moc/cell.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/surface.hpp>
#include <moc/pin_cell.hpp>
#include <moc/direction.hpp>
#include <moc/vector.hpp>
#include <transport_xs.hpp>

#include <xtensor/xarray.hpp>

#include <memory>
#include <optional>
#include <variant>
#include <vector>

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

    void append_fsrs(std::vector<FlatSourceRegion*>& fsrs);

    double trace_segments(Vector& r, const Direction& u,
                          std::vector<Segment>& segments);
    
    Tile clone() const {
      Tile t_out;

      if (this->c2d) {
        t_out.c2d = this->c2d->clone();
      }

      if (this->cell) {
        t_out.cell = this->cell->clone();
      }

      return t_out;
    }
  };

  struct TileIndex {
    std::size_t i, j;
  };

  Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
              const std::vector<std::shared_ptr<Surface>>& y_bounds);

  Cartesian2D(const std::vector<double>& dx, const std::vector<double>& dy);

  std::shared_ptr<Cartesian2D> clone() const;

  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }

  double dx() const { return x_bounds_.back()->x0() - x_bounds_.front()->x0(); }
  double dy() const { return y_bounds_.back()->y0() - y_bounds_.front()->y0(); }

  std::optional<TileIndex> get_tile_index(const Vector& r,
                                          const Direction& u) const;

  Vector get_tile_center(const TileIndex& ti) const;

  const Tile& tile(const TileIndex& ti) const;
  void set_tile(const TileIndex& ti, const std::shared_ptr<Cartesian2D>& c2d);
  void set_tile(const TileIndex& ti, const std::shared_ptr<Cell>& cell);
  void set_tiles(const std::vector<TileFill>& fills);

  bool tiles_valid() const;

  std::shared_ptr<TransportXS> get_xs(const Vector& r,
                                      const Direction& u) const;

  FlatSourceRegion& get_fsr(const Vector& r, const Direction& u);
  const FlatSourceRegion& get_fsr(const Vector& r, const Direction& u) const;

  std::size_t num_fsrs() const;

  void append_fsrs(std::vector<FlatSourceRegion*>& fsrs);

  double trace_segments(Vector& r, const Direction& u,
                        std::vector<Segment>& segments);

  double x_min() const { return x_bounds_.front()->x0(); }
  double x_max() const { return x_bounds_.back()->x0(); }

  double y_min() const { return y_bounds_.front()->y0(); }
  double y_max() const { return y_bounds_.back()->y0(); }

 private:
  std::vector<std::shared_ptr<Surface>> x_bounds_;
  std::vector<std::shared_ptr<Surface>> y_bounds_;
  xt::xarray<Tile> tiles_;
  std::size_t nx_, ny_;

  Cartesian2D() = default;

  void check_tile_index(const TileIndex& ti) const;
  std::pair<double, double> tile_dx_dy(const TileIndex& ti) const;
};

}  // namespace scarabee

#endif
