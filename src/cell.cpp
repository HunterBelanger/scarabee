#include <moc/cell.hpp>
#include <utils/scarabee_exception.hpp>

Cell::Cell(double xmin, double xmax, double ymin, double ymax)
    : fsrs_(), xss_(), x_min_(xmin), x_max_(xmax), y_min_(ymin), y_max_(ymax) {
  check_values();
}

void Cell::check_values() const {
  if (x_min_ >= x_max_) {
    throw ScarabeeException("x_min must be < x_max.");
  }

  if (y_min_ >= y_max_) {
    throw ScarabeeException("y_min must be < y_max.");
  }
}
