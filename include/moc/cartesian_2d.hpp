#ifndef CARTESIAN_2D_H
#define CARTESIAN_2D_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <transport_xs.hpp>

#include "cell.hpp"
#include "surface.hpp"
#include "vector.hpp"

#include <memory>
#include <vector>

class Cartesian2D;

class Cartesian2D {
 public:
  struct Tile {
    std::shared_ptr<Cartesian2D> c2d;
    std::shared_ptr<Cell> cell;
  };

  Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
              const std::vector<std::shared_ptr<Surface>>& y_bounds);

  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }

  const Tile& tile(std::size_t i, std::size_t j) const;

  void set_tile_PinCell(std::size_t i, std::size_t j,
                        const Vector& pin_origin,
                        const std::vector<double>& rads,
                        const std::vector<std::shared_ptr<TransportXS>>& mats);

 private:
  std::vector<std::shared_ptr<Surface>> x_bounds_;
  std::vector<std::shared_ptr<Surface>> y_bounds_;
  std::vector<Tile> tiles_;
  std::size_t nx_, ny_;
};

#endif