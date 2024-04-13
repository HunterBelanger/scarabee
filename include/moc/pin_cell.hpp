#ifndef PIN_CELL_H
#define PIN_CELL_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <transport_xs.hpp>

#include <memory>

namespace scarabee {

class PinCell : public Cell {
 public:
  PinCell(const std::vector<double>& rads,
          const std::vector<std::shared_ptr<TransportXS>>& mats,
          double dx, double dy);

 private:
  std::vector<double> mat_radii_;
  std::vector<std::shared_ptr<TransportXS>> mats_;
  std::vector<std::shared_ptr<Surface>> radii_;
  std::shared_ptr<Surface> xm_, pd_, ym_, nd_;


  void build();
};

}  // namespace scarabee

#endif
