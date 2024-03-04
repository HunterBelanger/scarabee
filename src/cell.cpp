#include <moc/cell.hpp>
#include <utils/scarabee_exception.hpp>

#include <cmath>
#include <sstream>

Cell::Cell(std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
           std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax)
    : fsrs_(), x_min_(xmin), x_max_(xmax), y_min_(ymin), y_max_(ymax) {}

std::vector<Segment> Cell::trace_segments(Vector& r, const Direction& u) {
  std::vector<Segment> segments;

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

  return segments;
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
  if (this->inside(r, u) == false) {
    throw ScarabeeException("Must be inside Cell to request FSR.");
  }

  for (auto& fsr : fsrs_) {
    if (fsr.inside(r, u)) return fsr;
  }

  std::stringstream mssg;
  mssg << "Could not find FSR at r = " << r << ", u = " << u << ".";
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return fsrs_.front();
}

const FlatSourceRegion& Cell::get_fsr(const Vector& r,
                                      const Direction& u) const {
  if (this->inside(r, u) == false) {
    throw ScarabeeException("Must be inside Cell to request FSR.");
  }

  for (const auto& fsr : fsrs_) {
    if (fsr.inside(r, u)) return fsr;
  }

  std::stringstream mssg;
  mssg << "Could not find FSR at r = " << r << ", u = " << u << ".";
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return fsrs_.front();
}
