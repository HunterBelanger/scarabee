#include <moc/simple_pin_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <algorithm>

namespace scarabee {

SimplePinCell::SimplePinCell(const std::vector<double>& mat_rads,
                 const std::vector<std::shared_ptr<TransportXS>>& mats,
                 double dx, double dy)
    : Cell(dx, dy),
      mat_radii_(mat_rads),
      mats_(mats),
      radii_() {
  this->build();
}

std::shared_ptr<Cell> SimplePinCell::clone() const {
  return std::make_shared<SimplePinCell>(*this);
}

void SimplePinCell::build() {
  // Clear vectors
  radii_.clear();
  fsrs_.clear();

  // Get variables
  const double x0_ = 0.;
  const double y0_ = 0.;
  const double dx = (x_max_->x0() - x_min_->x0());
  const double dy = (y_max_->y0() - y_min_->y0());

  // Make sure the numbers for mats and rads is coherent
  if (mat_radii_.size() + 1 != mats_.size()) {
    auto mssg = "Must have one more material than material radii.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure we have at least 2 mats and one radii
  if (mat_radii_.size() == 0) {
    auto mssg = "Must have at least one radiius.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure radii are sorted
  if (std::is_sorted(mat_radii_.begin(), mat_radii_.end()) == false) {
    auto mssg = "All radii must be sorted.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure radii all > 0
  if (mat_radii_[0] <= 0.) {
    auto mssg = "All radii must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure mats aren't nullptr
  for (const auto& mat : mats_) {
    if (!mat) {
      auto mssg = "Found TransportXS which was nullptr.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Make sure mats all have same ngroups
  std::size_t ngroups = mats_.front()->ngroups();
  for (std::size_t mi = 0; mi < mats_.size(); mi++) {
    if (mats_[mi]->ngroups() != ngroups) {
      auto mssg = "Materials have different numbers of groups.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Make sure largest radius doesn't intersect the cell walls
  if (x0_ - mat_radii_.back() < x_min_->x0() ||
      x0_ + mat_radii_.back() > x_max_->x0() ||
      y0_ - mat_radii_.back() < y_min_->y0() ||
      y0_ + mat_radii_.back() > y_max_->y0()) {
    auto mssg = "Outer material radius may not intersect cell wall.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, make all annular regions
  for (std::size_t i = 0; i < mat_radii_.size(); i++) {
    double r = mat_radii_[i];
    radii_.push_back(std::make_shared<Surface>());
    radii_.back()->type() = Surface::Type::Cylinder;
    radii_.back()->x0() = x0_;
    radii_.back()->y0() = y0_;
    radii_.back()->r() = r;

    double vol = PI * r * r;
    if (i > 0) {
      // Calculate volume of a ring.
      vol -= PI * mat_radii_[i - 1] * mat_radii_[i - 1];
    }

    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }
  }

  // We now make the outer region 
  fsrs_.emplace_back();
  fsrs_.back().volume() = dx*dy - (PI * mat_radii_.back() * mat_radii_.back());
  fsrs_.back().xs() = mats_.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});
}

}  // namespace scarabee
