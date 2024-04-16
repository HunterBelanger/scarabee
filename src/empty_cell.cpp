#include <moc/empty_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <algorithm>

namespace scarabee {

EmptyCell::EmptyCell(const std::shared_ptr<TransportXS>& mat, double dx, double dy)
    : Cell(dx, dy),
      mat_(mat) {
  // Just a single FSR
  fsrs_.emplace_back();
  fsrs_.back().volume() = dx*dy;
  fsrs_.back().xs() = mat_;
  fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});
}

std::shared_ptr<Cell> EmptyCell::clone() const {
  return std::make_shared<EmptyCell>(*this);
}

}  // namespace scarabee
