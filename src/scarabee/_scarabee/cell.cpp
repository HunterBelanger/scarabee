#include <moc/cell.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <sstream>

namespace scarabee {

Cell::Cell(double dx, double dy)
    : fsrs_(),
      x_min_(nullptr),
      y_min_(nullptr),
      x_max_(nullptr),
      y_max_(nullptr) {
  // Check delta's
  if (dx <= 0.) {
    auto mssg = "Cell dx must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (dy <= 0.) {
    auto mssg = "Cell dy must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Build surfaces
  x_min_ = std::make_shared<XPlane>(-0.5 * dx);
  x_max_ = std::make_shared<XPlane>(0.5 * dx);

  y_min_ = std::make_shared<YPlane>(-0.5 * dy);
  y_max_ = std::make_shared<YPlane>(0.5 * dy);

  this->check_surfaces();
}

void Cell::check_surfaces() const {
  // Make sure we don't have a nullptr surface
  if (!x_min_) {
    auto mssg = "No x_min surface provided to Cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (!x_max_) {
    auto mssg = "No x_max surface provided to Cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (!y_min_) {
    auto mssg = "No y_min surface provided to Cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (!y_max_) {
    auto mssg = "No y_max surface provided to Cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure the surfaces are ordered
  if (x_min_->x0() >= x_max_->x0()) {
    auto mssg = "Cell surface x_min must be < x_max.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (y_min_->y0() >= y_max_->y0()) {
    auto mssg = "Cell surface y_min must be < y_max.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

std::vector<UniqueFSR> Cell::get_all_fsr_in_cell(const Vector& /*r*/,
                                                 const Direction& /*u*/) const {
  std::vector<UniqueFSR> out;
  out.reserve(fsrs_.size());

  for (const auto& fsr : fsrs_) {
    out.push_back({&fsr, 0});
  }

  return out;
}

std::set<std::size_t> Cell::get_all_fsr_ids() const {
  std::set<std::size_t> ids;

  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    ids.insert(fsrs_[i].id());
  }

  return ids;
}

std::size_t Cell::get_num_fsr_instances(std::size_t id) const {
  // Each cell should only have up to 1 instance of a FSR
  for (const auto& fsr : fsrs_) {
    if (fsr.id() == id) return 1;
  }

  return 0;
}

void Cell::fill_fsrs(
    std::map<std::size_t, const FlatSourceRegion*>& fsrs) const {
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    if (fsrs.find(fsrs_[i].id()) == fsrs.end()) {
      fsrs[fsrs_[i].id()] = &fsrs_[i];
    }
  }
}

}  // namespace scarabee
