#ifndef EMPTY_CELL_H
#define EMPTY_CELL_H

#include <moc/cell.hpp>
#include <moc/surface.hpp>
#include <data/cross_section.hpp>

#include <memory>

namespace scarabee {

class EmptyCell : public Cell {
 public:
  EmptyCell(const std::shared_ptr<CrossSection>& mat, double dx, double dy);

 private:
  std::shared_ptr<CrossSection> mat_;
};

}  // namespace scarabee

#endif
