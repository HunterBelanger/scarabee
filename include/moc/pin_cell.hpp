#ifndef PIN_CELL_H
#define PIN_CELL_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <transport_xs.hpp>

#include <memory>

class PinCell : public Cell {
 public:
  PinCell(const std::vector<double>& rads,
          const std::vector<std::shared_ptr<TransportXS>>& mats,
          std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
          std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax);

  void build();

 private:
  std::vector<double> mat_radii_;
  std::vector<std::shared_ptr<TransportXS>> mats_;
  std::vector<std::shared_ptr<Surface>> radii_;
  std::shared_ptr<Surface> xm_, pd_, ym_, nd_;
  double x0_, y0_; // Origin of cell and rings
};

#endif