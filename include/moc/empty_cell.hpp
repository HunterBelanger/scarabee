#ifndef EMPTY_CELL_H
#define EMPTY_CELL_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <transport_xs.hpp>

#include <memory>

namespace scarabee {

class EmptyCell : public Cell {
 public:
  EmptyCell(const std::shared_ptr<TransportXS>& mat, double dx, double dy);

  std::shared_ptr<Cell> clone() const override final;

 private:
  std::shared_ptr<TransportXS> mat_;
};

}  // namespace scarabee

#endif
