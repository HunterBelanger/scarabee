#include <moc/cell.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <sstream>

namespace scarabee {

Cell::Cell(double dx, double dy)
    : fsrs_(),
      x_min_(nullptr),
      x_max_(nullptr),
      y_min_(nullptr),
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
  x_min_ = std::make_shared<Surface>();
  x_min_->type() = Surface::Type::XPlane;
  x_min_->x0() = -0.5 * dx;

  x_max_ = std::make_shared<Surface>();
  x_max_->type() = Surface::Type::XPlane;
  x_max_->x0() = 0.5 * dx;

  y_min_ = std::make_shared<Surface>();
  y_min_->type() = Surface::Type::YPlane;
  y_min_->y0() = -0.5 * dy;

  y_max_ = std::make_shared<Surface>();
  y_max_->type() = Surface::Type::YPlane;
  y_max_->y0() = 0.5 * dy;

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

bool Cell::inside(const Vector& r, const Direction& u) const {
  if (x_min_->side(r, u) == Surface::Side::Negative) return false;
  if (y_min_->side(r, u) == Surface::Side::Negative) return false;
  if (x_max_->side(r, u) == Surface::Side::Positive) return false;
  if (y_max_->side(r, u) == Surface::Side::Positive) return false;
  return true;
}

double Cell::distance(const Vector& r, const Direction& u) const {
  const double x_min_dist = x_min_->distance(r, u);
  const double x_max_dist = x_max_->distance(r, u);
  const double y_min_dist = y_min_->distance(r, u);
  const double y_max_dist = y_max_->distance(r, u);
  return std::min(std::min(x_min_dist, x_max_dist),
                  std::min(y_min_dist, y_max_dist));
}

UniqueFSR Cell::get_fsr(const Vector& r, const Direction& u) const {
  std::stringstream mssg;
  if (this->inside(r, u) == false) {
    mssg << "Could not find FSR at r = " << r << ", u = " << u << ".\n";
    mssg << "Position r and direction u are not inside the cell.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  for (const auto& fsr : fsrs_) {
    if (fsr.inside(r, u)) return {&fsr, 0};
  }

  mssg << "Could not find FSR at r = " << r << ", u = " << u << ".";
  spdlog::error(mssg.str());
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return {&fsrs_.front(), 0};
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
    if (fsr.id() == id)
      return 1;
  }

  return 0;
}

}  // namespace scarabee
