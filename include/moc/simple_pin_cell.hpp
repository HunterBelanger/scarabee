#ifndef SIMPLE_PIN_CELL_H
#define SIMPLE_PIN_CELL_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <transport_xs.hpp>

#include <memory>

namespace scarabee {

class SimplePinCell : public Cell {
 public:
  SimplePinCell(const std::vector<double>& rads, const std::vector<std::shared_ptr<TransportXS>>& mats, double dx,
          double dy);

 private:
  std::vector<double> mat_radii_;
  std::vector<std::shared_ptr<TransportXS>> mats_;
  std::vector<std::shared_ptr<Surface>> radii_;

  void build();
};

}  // namespace scarabee

#endif
