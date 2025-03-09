#include <moc/bwr_corner_pin_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <algorithm>
#include <cmath>
#include <map>

namespace scarabee {

BWRCornerPinCell::BWRCornerPinCell(
    const std::vector<double>& pin_rads,
    const std::vector<std::shared_ptr<CrossSection>>& pin_mats,
    std::shared_ptr<CrossSection> inner_mod, double inner_gap, double box_width,
    std::shared_ptr<CrossSection> box_mat,
    std::shared_ptr<CrossSection> outer_mod, double rc, double dx, double dy,
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
      xm_(nullptr),
      pd_(nullptr),
      ym_(nullptr),
      nd_(nullptr),
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

void BWRCornerPinCell::build_I() {
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

  // Create the FSRs for the outer moderator
  // position of outer box corner
  const Vector r_obc(Rcx + (rc_ + box_width_) / SQRT_2,
                     Rcy + (rc_ + box_width_) / SQRT_2);
  std::shared_ptr<Surface> robc_xp = std::make_shared<XPlane>(r_obc.x());
  std::shared_ptr<Surface> robc_yp = std::make_shared<YPlane>(r_obc.y());

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        robc_xp->integrate_y(r_obc.y(), ymax, Surface::Side::Positive) -
        box_outer->integrate_y(r_obc.y(), ymax, Surface::Side::Positive);
    fsr.xs() = outer_mod_;
    fsr.tokens() = {{box_outer, Surface::Side::Positive},
                    {robc_xp, Surface::Side::Negative},
                    {y_max_, Surface::Side::Negative}};
  }

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = outer_mod_;
    fsr.volume() = 0.5 * (xmax - r_obc.x()) * (ymax - r_obc.y());
    fsr.tokens() = {{robc_xp, Surface::Side::Positive},
                    {pd_, Surface::Side::Positive},
                    {y_max_, Surface::Side::Negative}};
  }

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = outer_mod_;
    fsr.volume() = 0.5 * (xmax - r_obc.x()) * (ymax - r_obc.y());
    fsr.tokens() = {{robc_yp, Surface::Side::Positive},
                    {pd_, Surface::Side::Negative},
                    {x_max_, Surface::Side::Negative}};
  }

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.volume() =
        robc_yp->integrate_x(r_obc.x(), xmax, Surface::Side::Positive) -
        box_outer->integrate_x(r_obc.x(), xmax, Surface::Side::Positive);
    fsr.xs() = outer_mod_;
    fsr.tokens() = {{box_outer, Surface::Side::Positive},
                    {robc_yp, Surface::Side::Negative},
                    {x_max_, Surface::Side::Negative}};
  }

  // Create the FSRs for the box
  std::shared_ptr<Surface> box_xp = std::make_shared<XPlane>(Rcx);
  std::shared_ptr<Surface> box_yp = std::make_shared<YPlane>(Rcy);

  if (rc_ < dx) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = box_mat_;
    fsr.volume() = box_width_ * (Rcx - xmin);
    fsr.tokens() = {{box_outer, Surface::Side::Negative},
                    {box_inner, Surface::Side::Positive},
                    {box_xp, Surface::Side::Negative}};
  }

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = box_mat_;
    fsr.volume() = 0.5 * 0.25 * PI * (rc_ + box_width_) *
                   (rc_ + box_width_);            // Larger cylinder
    fsr.volume() -= 0.5 * 0.25 * PI * rc_ * rc_;  // Subtract inner cylinder
    fsr.tokens() = {{box_outer, Surface::Side::Negative},
                    {box_inner, Surface::Side::Positive},
                    {box_xp, Surface::Side::Positive},
                    {pd_, Surface::Side::Positive}};
  }

  {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = box_mat_;
    fsr.volume() = 0.5 * 0.25 * PI * (rc_ + box_width_) *
                   (rc_ + box_width_);            // Larger cylinder
    fsr.volume() -= 0.5 * 0.25 * PI * rc_ * rc_;  // Subtract inner cylinder
    fsr.tokens() = {{box_outer, Surface::Side::Negative},
                    {box_inner, Surface::Side::Positive},
                    {pd_, Surface::Side::Negative},
                    {box_yp, Surface::Side::Positive}};
  }

  if (rc_ < dy) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = box_mat_;
    fsr.volume() = box_width_ * (Rcy - ymin);
    fsr.tokens() = {{box_outer, Surface::Side::Negative},
                    {box_inner, Surface::Side::Positive},
                    {box_yp, Surface::Side::Negative}};
  }

  // Create FSRs for the inner moderator gaps
  std::shared_ptr<Surface> gap_yp =
      std::make_shared<YPlane>(ymax - box_width_ - inner_gap_);
  std::shared_ptr<Surface> gap_xp =
      std::make_shared<XPlane>(xmax - box_width_ - inner_gap_);

  if (inner_gap_ > 0. && rc_ < dx) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = inner_mod_;
    fsr.volume() = inner_gap_ * (Rcx - xmin);
    fsr.tokens() = {{box_inner, Surface::Side::Negative},
                    {gap_yp, Surface::Side::Positive},
                    {box_xp, Surface::Side::Negative}};
  }

  if (inner_gap_ > 0.) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = inner_mod_;
    fsr.volume() =
        box_inner->integrate_y(ymax - box_width_ - inner_gap_,
                               ymax - box_width_, Surface::Side::Positive) -
        box_xp->integrate_y(ymax - box_width_ - inner_gap_, ymax - box_width_,
                            Surface::Side::Negative);
    fsr.tokens() = {{box_inner, Surface::Side::Negative},
                    {gap_yp, Surface::Side::Positive},
                    {box_xp, Surface::Side::Positive}};
  }

  if (inner_gap_ > 0.) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = inner_mod_;
    fsr.volume() =
        box_inner->integrate_x(xmax - box_width_ - inner_gap_,
                               xmax - box_width_, Surface::Side::Positive) -
        box_yp->integrate_x(xmax - box_width_ - inner_gap_, xmax - box_width_,
                            Surface::Side::Negative);
    fsr.tokens() = {{box_inner, Surface::Side::Negative},
                    {gap_xp, Surface::Side::Positive},
                    {box_yp, Surface::Side::Positive}};
  }

  if (inner_gap_ > 0. && rc_ < dy) {
    fsrs_.emplace_back();
    auto& fsr = fsrs_.back();
    fsr.xs() = inner_mod_;
    fsr.volume() = inner_gap_ * (Rcy - ymin);
    fsr.tokens() = {{box_inner, Surface::Side::Negative},
                    {gap_xp, Surface::Side::Positive},
                    {box_yp, Surface::Side::Negative}};
  }

  // Create the FSRs for the inner moderator pin-cell region
  // First get the radius of the moderator around pin
  const double pin_mod_rad = Rpx - xmin;
  surfs_.push_back(std::make_shared<Cylinder>(Rpx, Rpy, pin_mod_rad));

  // We need to check if the box will intersect the moderator ring around
  // the fuel pin. If so, we add an extra surface.
  const double Rx = 0.5 * (dx - box_width_ - inner_gap_) + inner_gap_;
  const double Ry = 0.5 * (dy - box_width_ - inner_gap_) + inner_gap_;
  const double R = std::sqrt(Rx * Rx + Ry * Ry);
  const double corner_test_criteria = (R - pin_mod_rad) / (std::sqrt(2.) - 1.);

  // Make FSRs for moderator ring around fuel pin
  {
    const double vol =
        0.125 * PI *
        (pin_mod_rad * pin_mod_rad - pin_radii_.back() * pin_radii_.back());
    const std::size_t i = surfs_.size() - 1;
    
    // NNE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {xm_, Surface::Side::Positive},
                      {pd_, Surface::Side::Positive}};
      if (rc_ > corner_test_criteria) {
        fsr.tokens().push_back({box_inner, Surface::Side::Negative});
        fsr.volume() = 0.;
      }
    }

    // NEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {ym_, Surface::Side::Positive},
                      {pd_, Surface::Side::Negative}};
      if (rc_ > corner_test_criteria) {
        fsr.tokens().push_back({box_inner, Surface::Side::Negative});
        fsr.volume() = 0.;
      }
    }

    // SEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {ym_, Surface::Side::Negative},
                      {nd_, Surface::Side::Positive}};
      if (rc_ > corner_test_criteria) {
        fsr.tokens().push_back({box_inner, Surface::Side::Negative});
        fsr.volume() = 0.;
      }
    }

    // SSE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {xm_, Surface::Side::Positive},
                      {nd_, Surface::Side::Negative}};
    }

    // SSW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {xm_, Surface::Side::Negative},
                      {pd_, Surface::Side::Negative}};
    }

    // SWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {ym_, Surface::Side::Negative},
                      {pd_, Surface::Side::Positive}};
    }

    // NWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {ym_, Surface::Side::Positive},
                      {nd_, Surface::Side::Negative}};
    }

    // NNW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_[i], Surface::Side::Negative},
                      {surfs_[i - 1], Surface::Side::Positive},
                      {xm_, Surface::Side::Negative},
                      {nd_, Surface::Side::Positive}};
      if (rc_ > corner_test_criteria) {
        fsr.tokens().push_back({box_inner, Surface::Side::Negative});
        fsr.volume() = 0.;
      }
    }
  }

  // Make the last FSRs for moderator outside the moderator ring
  {
    // NNE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {box_inner, Surface::Side::Negative},
                      {gap_yp, Surface::Side::Negative},
                      {xm_, Surface::Side::Positive},
                      {pd_, Surface::Side::Positive}};
    }

    // NEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {box_inner, Surface::Side::Negative},
                      {gap_xp, Surface::Side::Negative},
                      {ym_, Surface::Side::Positive},
                      {pd_, Surface::Side::Negative}};
    }

    // SEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {box_inner, Surface::Side::Negative},
                      {gap_xp, Surface::Side::Negative},
                      {ym_, Surface::Side::Negative},
                      {nd_, Surface::Side::Positive}};
    }

    // SSE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {y_min_, Surface::Side::Positive},
                      {xm_, Surface::Side::Positive},
                      {nd_, Surface::Side::Negative}};
    }

    // SSW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {y_min_, Surface::Side::Positive},
                      {xm_, Surface::Side::Negative},
                      {pd_, Surface::Side::Negative}};
    }

    // SWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {x_min_, Surface::Side::Positive},
                      {ym_, Surface::Side::Negative},
                      {pd_, Surface::Side::Positive}};
    }

    // NWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {x_min_, Surface::Side::Positive},
                      {ym_, Surface::Side::Positive},
                      {nd_, Surface::Side::Negative}};
    }

    // NNW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = 0.;
      fsr.xs() = inner_mod_;
      fsr.tokens() = {{surfs_.back(), Surface::Side::Positive},
                      {box_inner, Surface::Side::Negative},
                      {gap_yp, Surface::Side::Negative},
                      {xm_, Surface::Side::Negative},
                      {nd_, Surface::Side::Positive}};
    }
  }
  
  this->estimate_unknown_areas_cull_fsrs();
}

void BWRCornerPinCell::build_II() {
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

void BWRCornerPinCell::build_III() {
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

void BWRCornerPinCell::build_IV() {
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

void BWRCornerPinCell::build_pin(const double Rpx, const double Rpy) {
  // Make the plane surfaces for the angular divisions
  xm_ = std::make_shared<XPlane>(Rpx);
  ym_ = std::make_shared<YPlane>(Rpy);
  pd_ = std::make_shared<Plane>(-1., 1., Rpy - Rpx);
  nd_ = std::make_shared<Plane>(1., 1., Rpy + Rpx);

  // Make all annular regions
  for (std::size_t i = 0; i < pin_radii_.size(); i++) {
    const double r = pin_radii_[i];
    surfs_.push_back(std::make_shared<Cylinder>(Rpx, Rpy, r));

    // Make 8 FSRs
    double vol = PI * r * r;
    if (i > 0) {
      // Calculate volume of a ring.
      vol -= PI * pin_radii_[i - 1] * pin_radii_[i - 1];
    }
    // Get 1/8th the volume of the full ring, for each angular segment.
    vol /= 8.;

    // NNE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({xm_, Surface::Side::Positive});
      fsr.tokens().push_back({pd_, Surface::Side::Positive});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // NEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({ym_, Surface::Side::Positive});
      fsr.tokens().push_back({pd_, Surface::Side::Negative});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // SEE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({ym_, Surface::Side::Negative});
      fsr.tokens().push_back({nd_, Surface::Side::Positive});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // SSE
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({xm_, Surface::Side::Positive});
      fsr.tokens().push_back({nd_, Surface::Side::Negative});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // SSW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({xm_, Surface::Side::Negative});
      fsr.tokens().push_back({pd_, Surface::Side::Negative});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // SWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({ym_, Surface::Side::Negative});
      fsr.tokens().push_back({pd_, Surface::Side::Positive});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // NWW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({ym_, Surface::Side::Positive});
      fsr.tokens().push_back({nd_, Surface::Side::Negative});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }

    // NNW
    {
      fsrs_.emplace_back();
      auto& fsr = fsrs_.back();
      fsr.volume() = vol;
      fsr.xs() = pin_mats_[i];
      fsr.tokens().push_back({surfs_.back(), Surface::Side::Negative});
      fsr.tokens().push_back({xm_, Surface::Side::Negative});
      fsr.tokens().push_back({nd_, Surface::Side::Positive});
      if (i > 0) {
        fsr.tokens().push_back({surfs_[i - 1], Surface::Side::Positive});
      }
    }
  }
}

void BWRCornerPinCell::estimate_unknown_areas_cull_fsrs() {
  const double ymin = y_min_->y0();
  const double xmin = x_min_->x0();
  const double xmax = x_max_->x0();
  const double dx = xmax - xmin;
  std::size_t ntracks = static_cast<std::size_t>(std::round(dx / 0.0001));
  const double dt = dx / static_cast<double>(ntracks);

  // Initialize empty areas
  std::map<std::size_t, double> areas;
  for (auto& fsr : fsrs_) {
    if (fsr.volume() == 0.) {
      areas[fsr.id()] = 0.;
    }
  }

  // Do track integration
  const Direction u(0., 1.);
  for (std::size_t t = 0; t < ntracks; t++) {
    Vector r(xmin + (0.5 + static_cast<double>(t))*dt, ymin);

    auto ufsr = this->get_fsr(r, u);
    while (ufsr.fsr) {
      // Get length
      const double d = ufsr.fsr->distance(r, u);

      // Tally to area
      const auto fsr_id = ufsr.fsr->id();
      auto areas_iter = areas.find(fsr_id);
      if (areas_iter != areas.end()) {
        areas_iter->second += dt * d;
      }

      // Update position and ufsr
      r = Vector(r.x(), r.y()+d);
      ufsr = this->get_fsr(r, u);
    }
  }

  // Assign all areas
  for (auto& areas_iter : areas) {
    for (auto& fsr : fsrs_) {
      if (areas_iter.first == fsr.id()) {
        fsr.volume() = areas_iter.second;
        break;
      }
    }
  }

  // Remove FSRs that still have zero volume. Apparently they aren't visible !
  // Initially tried using std::remove_if, but this lead to odd artifacts at
  // the end of the fsrs_ vector. Not sure why, so I just made this.
  for (auto fsr_itr = fsrs_.begin(); fsr_itr != fsrs_.end(); fsr_itr++) {
    if (fsr_itr->volume() == 0.) {
      fsr_itr = fsrs_.erase(fsr_itr);
      fsr_itr--;
    }
  }
}

}  // namespace scarabee
