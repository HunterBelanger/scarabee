#include <moc/cartesian_2d.hpp>
#include <utils/scarabee_exception.hpp>

#include <algorithm>

Cartesian2D::Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
                         const std::vector<std::shared_ptr<Surface>>& y_bounds)
    : x_bounds_(x_bounds), y_bounds_(y_bounds), tiles_(), nx_(), ny_() {
  // First, make sure we have at least 2 bounds in each direction
  if (x_bounds_.size() < 2) {
    throw ScarabeeException("Must provide at least 2 x-bounds");
  } else if (y_bounds_.size() < 2) {
    throw ScarabeeException("Must provide at least 2 y-bounds");
  }

  // Make sure bounds are sorted
  if (std::is_sorted(x_bounds_.begin(), x_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->x0() < b->x0();
                     }) == false) {
    throw ScarabeeException("Provided x-bounds are not sorted.");
  }
  
  if (std::is_sorted(y_bounds_.begin(), y_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->y0() < b->y0();
                     }) == false) {
    throw ScarabeeException("Provided y-bounds are not sorted.");
  }

  // Get sizes and allocate tiles array
  nx_ = x_bounds_.size();
  ny_ = y_bounds_.size();

  // All tiles start as uninitialized
  tiles_.resize(nx_ * ny_, {nullptr, std::nullopt});
}

bool Cartesian2D::tiles_valid() const {
  for (const auto& t : tiles_) {
    if (t.c2d && t.cell) {
      return false;
    } else if (!t.c2d & !t.cell) {
      return false;
    } else if (t.c2d) {
      if (t.c2d->tiles_valid() == false) return false;
    } else {
      if (!t.cell) return false;
    }
  }

  return true;
}