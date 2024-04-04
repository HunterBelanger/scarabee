#include <moc/cell.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <cmath>
#include <sstream>

Cell::Cell(std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
           std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax)
    : fsrs_(), x_min_(xmin), x_max_(xmax), y_min_(ymin), y_max_(ymax) {
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

void Cell::set_x_min(const std::shared_ptr<Surface>& xm) {
  x_min_ = xm;
  this->check_surfaces();
}

void Cell::set_x_max(const std::shared_ptr<Surface>& xm) {
  x_max_ = xm;
  this->check_surfaces();
}

void Cell::set_y_min(const std::shared_ptr<Surface>& ym) {
  y_min_ = ym;
  this->check_surfaces();
}

void Cell::set_y_max(const std::shared_ptr<Surface>& ym) {
  y_max_ = ym;
  this->check_surfaces();
}

std::vector<Segment> Cell::trace_segments(Vector& r, const Direction& u) {
  std::vector<Segment> segments;
  this->trace_segments(r, u, segments);
  return segments;
}

void Cell::trace_segments(Vector& r, const Direction& u,
                          std::vector<Segment>& segments) {
  while (this->inside(r, u)) {
    // Get current FSR
    auto& fsr = this->get_fsr(r, u);

    // Get distance we can travel in FSR
    double distance = fsr.distance(r, u);

    // Add new segment
    segments.emplace_back(&fsr, distance);

    // Increment distance
    r = r + distance * u;
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

FlatSourceRegion& Cell::get_fsr(const Vector& r, const Direction& u) {
  std::stringstream mssg;
  if (this->inside(r, u) == false) {
    mssg << "Could not find FSR at r = " << r << ", u = " << u << ".\n";
    mssg << "Position r and direction u are not inside the cell.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  for (auto& fsr : fsrs_) {
    if (fsr.inside(r, u)) return fsr;
  }

  mssg << "Could not find FSR at r = " << r << ", u = " << u << ".";
  spdlog::error(mssg.str());
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return fsrs_.front();
}

const FlatSourceRegion& Cell::get_fsr(const Vector& r,
                                      const Direction& u) const {
  std::stringstream mssg;
  if (this->inside(r, u) == false) {
    mssg << "Could not find FSR at r = " << r << ", u = " << u << ".\n";
    mssg << "Position r and direction u are not inside the cell.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  for (const auto& fsr : fsrs_) {
    if (fsr.inside(r, u)) return fsr;
  }

  mssg << "Could not find FSR at r = " << r << ", u = " << u << ".";
  spdlog::error(mssg.str());
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return fsrs_.front();
}
