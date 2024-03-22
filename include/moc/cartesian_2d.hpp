#ifndef CARTESIAN_2D_H
#define CARTESIAN_2D_H

#include <moc/cell.hpp>
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

// Must pre-declare class so that we can have a pointer to it in the tile
class Cartesian2D;

using CellType = std::variant<PinCell>;

class Cartesian2D {
 public:
  struct Tile {
    std::shared_ptr<Cartesian2D> c2d;
    std::optional<CellType> cell;

    bool valid() const {
      if (c2d && cell) return false;
      return c2d || cell;
    }

    void trace_segments(Vector& r, const Direction& u,
                        std::vector<Segment>& segments);
  };

  struct TileIndex {
    std::size_t i, j;
  };

  Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
              const std::vector<std::shared_ptr<Surface>>& y_bounds);

  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }

  double x_min() const { return x_bounds_.front()->x0(); }
  double x_max() const { return x_bounds_.back()->x0(); }

  double y_min() const { return y_bounds_.front()->y0(); }
  double y_max() const { return y_bounds_.back()->y0(); }

  std::optional<TileIndex> get_tile_index(const Vector& r,
                                          const Direction& u) const;

  Tile& tile(const TileIndex& ti);
  const Tile& tile(const TileIndex& ti) const;

  bool tiles_valid() const;

  void trace_segments(Vector& r, const Direction& u,
                      std::vector<Segment>& segments);

 private:
  std::vector<std::shared_ptr<Surface>> x_bounds_;
  std::vector<std::shared_ptr<Surface>> y_bounds_;
  xt::xarray<Tile> tiles_;
  std::size_t nx_, ny_;
};

#endif
