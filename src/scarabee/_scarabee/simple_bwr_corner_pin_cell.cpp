#include <moc/simple_bwr_corner_pin_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <algorithm>
#include <cmath>

namespace scarabee {

SimpleBWRCornerPinCell::SimpleBWRCornerPinCell(
    const std::vector<double>& pin_rads,
    const std::vector<std::shared_ptr<CrossSection>>& pin_mats,
    double inner_gap, std::shared_ptr<CrossSection> inner_mod, double box_width,
    double rc, std::shared_ptr<CrossSection> box_mat,
    std::shared_ptr<CrossSection> outer_mod, double dx, double dy,
    BWRCornerType corner_type)
    : Cell(dx, dy),
      pin_radii_(pin_rads),
      pin_mats_(pin_mats),
      surfs_(),
      inner_mod_(inner_mod),
      inner_gap_(inner_gap),
      box_width_(box_width),
      box_mat_(box_mat),
      outer_mod_(outer_mod),
      rc_(rc),
      corner_type_(corner_type) {
  // Make sure we have at least one mat and one radii
  if (pin_radii_.size() == 0) {
    auto mssg = "Must have at least one radius.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (pin_mats_.size() != pin_radii_.size()) {
    auto mssg = "Number of pin radii and materials must agree.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure radii are sorted
  if (std::is_sorted(pin_radii_.begin(), pin_radii_.end()) == false) {
    auto mssg = "All radii must be sorted.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure radii all > 0
  if (pin_radii_[0] <= 0.) {
    auto mssg = "All radii must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure mats aren't nullptr
  for (const auto& mat : pin_mats_) {
    if (!mat) {
      auto mssg = "Found pin CrossSection which is None.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  if (!box_mat_) {
    auto mssg = "Box CrossSection is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (!inner_mod_) {
    auto mssg = "Inner moderator CrossSection is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (!outer_mod_) {
    auto mssg = "Outer moderator CrossSection is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure mats all have same ngroups
  std::size_t ngroups = pin_mats_.front()->ngroups();
  for (std::size_t mi = 0; mi < pin_mats_.size(); mi++) {
    if (pin_mats_[mi]->ngroups() != ngroups) {
      auto mssg = "Materials have different numbers of groups.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  if (inner_mod_->ngroups() != ngroups) {
    auto mssg = "Inner moderator has a different number of groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (box_mat_->ngroups() != ngroups) {
    auto mssg = "Box has a different number of groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (outer_mod_->ngroups() != ngroups) {
    auto mssg = "Outer moderator has a different number of groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure box has positive width
  if (box_width <= 0.) {
    auto mssg = "Box width must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure corner radius is >= 0.
  if (rc_ < 0.) {
    auto mssg = "Box corner radius >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure radius is smaller than the cell widths
  if (dx - (rc_ + box_width_) < 0.) {
    auto mssg = "Corner radius is too large for cell x-width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (dy - (rc_ + box_width_) < 0.) {
    auto mssg = "Corner radius is too large for cell y-width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure inner gap width isn't negative
  if (inner_gap_ < 0.) {
    auto mssg = "Inner gap width must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure inner gap width isn't too big !
  if (inner_gap_ > 0.) {
    // Computation performed in II corner frame
    // Position of upper left pin cell corner
    const double xpc = -0.5 * dx + box_width_ + inner_gap_;
    const double ypc = 0.5 * dy - box_width_ - inner_gap_;

    // Position of circle center
    const double Rcx = -0.5 * dx + rc_ + box_width_;
    const double Rcy = 0.5 * dy - (rc_ + box_width_);

    // Compute distance from circle center to the cell corner
    const double d =
        std::sqrt((Rcx - xpc) * (Rcx - xpc) + (Rcy - ypc) * (Rcy - ypc));
    if (d < rc_) {
      auto mssg = "Inner moderator gap is too wide for the cell.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Make sure the pins won't intersect the boundaries of the pin part
  if (pin_radii_.back() > 0.5 * (dx - box_width_ - inner_gap_)) {
    auto mssg = "Pin radius is too large for cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (pin_radii_.back() > 0.5 * (dy - box_width_ - inner_gap_)) {
    auto mssg = "Pin radius is too large for cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check if the pin will hit the rounded corner of the box.
  // Assume pin is at origin, and the "missing" corner would then be at
  const double Rx = 0.5 * (dx - box_width - inner_gap_) + inner_gap_;
  const double Ry = 0.5 * (dy - box_width - inner_gap_) + inner_gap_;
  const double R = std::sqrt(Rx * Rx + Ry * Ry);
  const double corner_test_criteria =
      (R - pin_radii_.back()) / (std::sqrt(2.) - 1.);
  if (rc_ > corner_test_criteria) {
    auto mssg =
        "Box corner radius is too large for given pin radius and cell size.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Now we can go ahead and build the cell
  switch (corner_type_) {
    case BWRCornerType::I:
      build_I();
      break;

    case BWRCornerType::II:
      build_II();
      break;

    case BWRCornerType::III:
      build_III();
      break;

    case BWRCornerType::IV:
      build_IV();
      break;
  }
}

void SimpleBWRCornerPinCell::build_I() {
  // The box corner is in the (+x,+y) corner of the cell.
  // Start be getting the position of the center of the pin.
  const double xmin = x_min_->x0();
  const double xmax = x_max_->x0();
  const double ymin = y_min_->y0();
  const double ymax = y_max_->y0();
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  const double pin_cell_dx = dx - box_width_ - inner_gap_;
  const double pin_cell_dy = dy - box_width_ - inner_gap_;

  // Compute pin center position in cell
  double Rpx = xmax - box_width_ - inner_gap_ - 0.5 * pin_cell_dx;
  double Rpy = ymax - box_width_ - inner_gap_ - 0.5 * pin_cell_dy;

  // Make the pin FSRs
  build_pin(Rpx, Rpy);

  // Compute center of the corner circle
  const double Rcx = xmax - rc_ - box_width_;
  const double Rcy = ymax - rc_ - box_width_;

  // Next, we make the two curved box surfaces.
  auto box_outer =
      std::make_shared<BWRCornerI>(xmin, xmax, ymin, ymax, rc_ + box_width_);
  auto box_inner = std::make_shared<BWRCornerI>(xmin, xmax - box_width_, ymin,
                                                ymax - box_width_, rc_);
  surfs_.push_back(box_outer);
  surfs_.push_back(box_inner);

  // Create the FSR for the outer moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() = (xmax - Rcx) * (ymax - Rcy);  // Box volume
    fsr.volume() -= 0.25 * PI * (rc_ + box_width_) *
                    (rc_ + box_width_);  // Subtract quarter circle
    fsr.xs() = outer_mod_;
    fsr.tokens().push_back({x_min_, Surface::Side::Positive});
    fsr.tokens().push_back({x_max_, Surface::Side::Negative});
    fsr.tokens().push_back({y_min_, Surface::Side::Positive});
    fsr.tokens().push_back({y_max_, Surface::Side::Negative});
    fsr.tokens().push_back({box_outer, Surface::Side::Positive});
  }

  // Create the FSR for the box
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        0.25 * PI * (rc_ + box_width_) * (rc_ + box_width_);  // Larger cylinder
    fsr.volume() -= 0.25 * PI * rc_ * rc_;      // Subtract inner cylinder
    fsr.volume() += box_width_ * (Rcy - ymin);  // Add first square
    fsr.volume() += box_width_ * (Rcx - xmin);  // Add second square
    fsr.xs() = box_mat_;
    fsr.tokens().push_back({box_outer, Surface::Side::Negative});
    fsr.tokens().push_back({box_inner, Surface::Side::Positive});
  }

  // Create the FSR for the inner moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    // Start with area of the rectangular portion.
    fsr.volume() = (dx - box_width_) * (dy - box_width_);
    // Remove the bit outside the cylinder
    fsr.volume() -= rc_ * rc_ * (1. - 0.25 * PI);
    // Remove area of pin
    if (pin_radii_.size() > 0)
      fsr.volume() -= PI * pin_radii_.back() * pin_radii_.back();
    fsr.xs() = inner_mod_;
    fsr.tokens().push_back({box_inner, Surface::Side::Negative});
    fsr.tokens().push_back(
        {surfs_[pin_radii_.size() - 1], Surface::Side::Positive});
  }
}

void SimpleBWRCornerPinCell::build_II() {
  // The box corner is in the (-x,+y) corner of the cell.
  // Start be getting the position of the center of the pin.
  const double xmin = x_min_->x0();
  const double xmax = x_max_->x0();
  const double ymin = y_min_->y0();
  const double ymax = y_max_->y0();
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  const double pin_cell_dx = dx - box_width_ - inner_gap_;
  const double pin_cell_dy = dy - box_width_ - inner_gap_;

  // Compute pin center position in cell
  double Rpx = xmin + box_width_ + inner_gap_ + 0.5 * pin_cell_dx;
  double Rpy = ymax - box_width_ - inner_gap_ - 0.5 * pin_cell_dy;

  // Make the pin FSRs
  build_pin(Rpx, Rpy);

  // Compute center of the corner circle
  const double Rcx = xmin + rc_ + box_width_;
  const double Rcy = ymax - rc_ - box_width_;

  // Next, we make the two curved box surfaces.
  auto box_outer =
      std::make_shared<BWRCornerII>(xmin, xmax, ymin, ymax, rc_ + box_width_);
  auto box_inner = std::make_shared<BWRCornerII>(xmin + box_width_, xmax, ymin,
                                                 ymax - box_width_, rc_);
  surfs_.push_back(box_outer);
  surfs_.push_back(box_inner);

  // Create the FSR for the outer moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() = (Rcx - xmin) * (ymax - Rcy);  // Box volume
    fsr.volume() -= 0.25 * PI * (rc_ + box_width_) *
                    (rc_ + box_width_);  // Subtract quarter circle
    fsr.xs() = outer_mod_;
    fsr.tokens().push_back({x_min_, Surface::Side::Positive});
    fsr.tokens().push_back({x_max_, Surface::Side::Negative});
    fsr.tokens().push_back({y_min_, Surface::Side::Positive});
    fsr.tokens().push_back({y_max_, Surface::Side::Negative});
    fsr.tokens().push_back({box_outer, Surface::Side::Positive});
  }

  // Create the FSR for the box
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        0.25 * PI * (rc_ + box_width_) * (rc_ + box_width_);  // Larger cylinder
    fsr.volume() -= 0.25 * PI * rc_ * rc_;      // Subtract inner cylinder
    fsr.volume() += box_width_ * (Rcy - ymin);  // Add first square
    fsr.volume() += box_width_ * (xmax - Rcx);  // Add second square
    fsr.xs() = box_mat_;
    fsr.tokens().push_back({box_outer, Surface::Side::Negative});
    fsr.tokens().push_back({box_inner, Surface::Side::Positive});
  }

  // Create the FSR for the inner moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    // Start with area of the rectangular portion.
    fsr.volume() = (dx - box_width_) * (dy - box_width_);
    // Remove the bit outside the cylinder
    fsr.volume() -= rc_ * rc_ * (1. - 0.25 * PI);
    // Remove area of pin
    if (pin_radii_.size() > 0)
      fsr.volume() -= PI * pin_radii_.back() * pin_radii_.back();
    fsr.xs() = inner_mod_;
    fsr.tokens().push_back({box_inner, Surface::Side::Negative});
    fsr.tokens().push_back(
        {surfs_[pin_radii_.size() - 1], Surface::Side::Positive});
  }
}

void SimpleBWRCornerPinCell::build_III() {
  // The box corner is in the (-x,-y) corner of the cell.
  // Start be getting the position of the center of the pin.
  const double xmin = x_min_->x0();
  const double xmax = x_max_->x0();
  const double ymin = y_min_->y0();
  const double ymax = y_max_->y0();
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  const double pin_cell_dx = dx - box_width_ - inner_gap_;
  const double pin_cell_dy = dy - box_width_ - inner_gap_;

  // Compute pin center position in cell
  double Rpx = xmin + box_width_ + inner_gap_ + 0.5 * pin_cell_dx;
  double Rpy = ymin + box_width_ + inner_gap_ + 0.5 * pin_cell_dy;

  // Make the pin FSRs
  build_pin(Rpx, Rpy);

  // Compute center of the corner circle
  const double Rcx = xmin + rc_ + box_width_;
  const double Rcy = ymin + rc_ + box_width_;

  // Next, we make the two curved box surfaces.
  auto box_outer =
      std::make_shared<BWRCornerIII>(xmin, xmax, ymin, ymax, rc_ + box_width_);
  auto box_inner = std::make_shared<BWRCornerIII>(xmin + box_width_, xmax,
                                                  ymin + box_width_, ymax, rc_);
  surfs_.push_back(box_outer);
  surfs_.push_back(box_inner);

  // Create the FSR for the outer moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() = (Rcx - xmin) * (Rcy - ymin);  // Box volume
    fsr.volume() -= 0.25 * PI * (rc_ + box_width_) *
                    (rc_ + box_width_);  // Subtract quarter circle
    fsr.xs() = outer_mod_;
    fsr.tokens().push_back({x_min_, Surface::Side::Positive});
    fsr.tokens().push_back({x_max_, Surface::Side::Negative});
    fsr.tokens().push_back({y_min_, Surface::Side::Positive});
    fsr.tokens().push_back({y_max_, Surface::Side::Negative});
    fsr.tokens().push_back({box_outer, Surface::Side::Positive});
  }

  // Create the FSR for the box
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        0.25 * PI * (rc_ + box_width_) * (rc_ + box_width_);  // Larger cylinder
    fsr.volume() -= 0.25 * PI * rc_ * rc_;      // Subtract inner cylinder
    fsr.volume() += box_width_ * (ymax - Rcy);  // Add first square
    fsr.volume() += box_width_ * (xmax - Rcx);  // Add second square
    fsr.xs() = box_mat_;
    fsr.tokens().push_back({box_outer, Surface::Side::Negative});
    fsr.tokens().push_back({box_inner, Surface::Side::Positive});
  }

  // Create the FSR for the inner moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    // Start with area of the rectangular portion.
    fsr.volume() = (dx - box_width_) * (dy - box_width_);
    // Remove the bit outside the cylinder
    fsr.volume() -= rc_ * rc_ * (1. - 0.25 * PI);
    // Remove area of pin
    if (pin_radii_.size() > 0)
      fsr.volume() -= PI * pin_radii_.back() * pin_radii_.back();
    fsr.xs() = inner_mod_;
    fsr.tokens().push_back({box_inner, Surface::Side::Negative});
    fsr.tokens().push_back(
        {surfs_[pin_radii_.size() - 1], Surface::Side::Positive});
  }
}

void SimpleBWRCornerPinCell::build_IV() {
  // The box corner is in the (+x,-y) corner of the cell.
  // Start be getting the position of the center of the pin.
  const double xmin = x_min_->x0();
  const double xmax = x_max_->x0();
  const double ymin = y_min_->y0();
  const double ymax = y_max_->y0();
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  const double pin_cell_dx = dx - box_width_ - inner_gap_;
  const double pin_cell_dy = dy - box_width_ - inner_gap_;

  // Compute pin center position in cell
  double Rpx = xmax - box_width_ - inner_gap_ - 0.5 * pin_cell_dx;
  double Rpy = ymin + box_width_ + inner_gap_ + 0.5 * pin_cell_dy;

  // Make the pin FSRs
  build_pin(Rpx, Rpy);

  // Compute center of the corner circle
  const double Rcx = xmax - rc_ - box_width_;
  const double Rcy = ymin + rc_ + box_width_;

  // Next, we make the two curved box surfaces.
  auto box_outer =
      std::make_shared<BWRCornerIV>(xmin, xmax, ymin, ymax, rc_ + box_width_);
  auto box_inner = std::make_shared<BWRCornerIV>(xmin, xmax - box_width_,
                                                 ymin + box_width_, ymax, rc_);
  surfs_.push_back(box_outer);
  surfs_.push_back(box_inner);

  // Create the FSR for the outer moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() = (xmax - Rcx) * (Rcy - ymin);  // Box volume
    fsr.volume() -= 0.25 * PI * (rc_ + box_width_) *
                    (rc_ + box_width_);  // Subtract quarter circle
    fsr.xs() = outer_mod_;
    fsr.tokens().push_back({x_min_, Surface::Side::Positive});
    fsr.tokens().push_back({x_max_, Surface::Side::Negative});
    fsr.tokens().push_back({y_min_, Surface::Side::Positive});
    fsr.tokens().push_back({y_max_, Surface::Side::Negative});
    fsr.tokens().push_back({box_outer, Surface::Side::Positive});
  }

  // Create the FSR for the box
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        0.25 * PI * (rc_ + box_width_) * (rc_ + box_width_);  // Larger cylinder
    fsr.volume() -= 0.25 * PI * rc_ * rc_;      // Subtract inner cylinder
    fsr.volume() += box_width_ * (ymax - Rcy);  // Add first square
    fsr.volume() += box_width_ * (Rcx - xmin);  // Add second square
    fsr.xs() = box_mat_;
    fsr.tokens().push_back({box_outer, Surface::Side::Negative});
    fsr.tokens().push_back({box_inner, Surface::Side::Positive});
  }

  // Create the FSR for the inner moderator
  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    // Start with area of the rectangular portion.
    fsr.volume() = (dx - box_width_) * (dy - box_width_);
    // Remove the bit outside the cylinder
    fsr.volume() -= rc_ * rc_ * (1. - 0.25 * PI);
    // Remove area of pin
    if (pin_radii_.size() > 0)
      fsr.volume() -= PI * pin_radii_.back() * pin_radii_.back();
    fsr.xs() = inner_mod_;
    fsr.tokens().push_back({box_inner, Surface::Side::Negative});
    fsr.tokens().push_back(
        {surfs_[pin_radii_.size() - 1], Surface::Side::Positive});
  }
}

void SimpleBWRCornerPinCell::build_pin(const double Rpx, const double Rpy) {
  // Make all annular regions
  for (std::size_t i = 0; i < pin_radii_.size(); i++) {
    const double r = pin_radii_[i];
    surfs_.push_back(std::make_shared<Cylinder>(Rpx, Rpy, r));

    double vol = PI * r * r;
    if (i > 0) {
      // Calculate volume of a ring.
      vol -= PI * pin_radii_[i - 1] * pin_radii_[i - 1];
    }

    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = pin_mats_[i];
    fsrs_.back().tokens().push_back({surfs_.back(), Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
    }
  }
}

}  // namespace scarabee
