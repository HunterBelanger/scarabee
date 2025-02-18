#include <moc/pin_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <algorithm>

namespace scarabee {

PinCell::PinCell(const std::vector<double>& mat_rads,
                 const std::vector<std::shared_ptr<CrossSection>>& mats,
                 double dx, double dy, PinCellType pin_type)
    : Cell(dx, dy),
      mat_radii_(mat_rads),
      mats_(mats),
      radii_(),
      xm_(),
      pd_(),
      ym_(),
      nd_(),
      pin_type_(pin_type) {
  this->build();
}

void PinCell::build() {
  // Clear vectors
  radii_.clear();
  fsrs_.clear();

  // Width of the cell in each direction
  const double dx = (x_max_->x0() - x_min_->x0());
  const double dy = (y_max_->y0() - y_min_->y0());

  // The center point of the pin
  double x0_ = 0.;
  double y0_ = 0.;

  if (pin_type_ == PinCellType::XN || pin_type_ == PinCellType::II ||
      pin_type_ == PinCellType::III) {
    x0_ = 0.5 * dx;
  } else if (pin_type_ == PinCellType::XP || pin_type_ == PinCellType::I ||
             pin_type_ == PinCellType::IV) {
    x0_ = -0.5 * dx;
  }

  if (pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III ||
      pin_type_ == PinCellType::IV) {
    y0_ = 0.5 * dy;
  } else if (pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I ||
             pin_type_ == PinCellType::II) {
    y0_ = -0.5 * dy;
  }

  // Create the 4 surfaces which give us our 8 angular sections
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::YN ||
      pin_type_ == PinCellType::YP) {
    xm_ = std::make_shared<XPlane>(x0_);
  } else if (pin_type_ == PinCellType::XP || pin_type_ == PinCellType::I ||
             pin_type_ == PinCellType::IV) {
    xm_ = x_min_;
  } else if (pin_type_ == PinCellType::XN || pin_type_ == PinCellType::II ||
             pin_type_ == PinCellType::III) {
    xm_ = x_max_;
  }

  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
      pin_type_ == PinCellType::XP) {
    ym_ = std::make_shared<YPlane>(y0_);
  } else if (pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I ||
             pin_type_ == PinCellType::II) {
    ym_ = y_min_;
  } else if (pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III ||
             pin_type_ == PinCellType::IV) {
    ym_ = y_max_;
  }

  pd_ = std::make_shared<Plane>(-1., 1., y0_ - x0_);

  nd_ = std::make_shared<Plane>(1., 1., y0_ + x0_);

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
      auto mssg = "Found CrossSection which was nullptr.";
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
  if (pin_type_ == PinCellType::Full) {
    if (x0_ - mat_radii_.back() < x_min_->x0() ||
        x0_ + mat_radii_.back() > x_max_->x0() ||
        y0_ - mat_radii_.back() < y_min_->y0() ||
        y0_ + mat_radii_.back() > y_max_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::XP) {
    if (x0_ + mat_radii_.back() > x_max_->x0() ||
        y0_ + mat_radii_.back() > y_max_->y0() ||
        y0_ - mat_radii_.back() < y_min_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::XN) {
    if (x0_ - mat_radii_.back() < x_min_->x0() ||
        y0_ + mat_radii_.back() > y_max_->y0() ||
        y0_ - mat_radii_.back() < y_min_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::YP) {
    if (y0_ + mat_radii_.back() > y_max_->y0() ||
        x0_ + mat_radii_.back() > x_max_->x0() ||
        x0_ - mat_radii_.back() < x_min_->x0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::YN) {
    if (y0_ - mat_radii_.back() < y_min_->y0() ||
        x0_ + mat_radii_.back() > x_max_->x0() ||
        x0_ - mat_radii_.back() < x_min_->x0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::I) {
    if (x0_ + mat_radii_.back() > x_max_->x0() ||
        y0_ + mat_radii_.back() > y_max_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::II) {
    if (x0_ - mat_radii_.back() < x_min_->x0() ||
        y0_ + mat_radii_.back() > y_max_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (pin_type_ == PinCellType::III) {
    if (x0_ - mat_radii_.back() < x_min_->x0() ||
        y0_ - mat_radii_.back() < y_min_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else {
    // pin_type_ == PinCellType::IV
    if (x0_ + mat_radii_.back() > x_max_->x0() ||
        y0_ - mat_radii_.back() < y_min_->y0()) {
      auto mssg = "Outer material radius may not intersect cell wall.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // First, make all annular regions
  for (std::size_t i = 0; i < mat_radii_.size(); i++) {
    double r = mat_radii_[i];
    radii_.push_back(std::make_shared<Cylinder>(x0_, y0_, r));

    // Make 8 FSRs
    double vol = PI * r * r;
    if (i > 0) {
      // Calculate volume of a ring.
      vol -= PI * mat_radii_[i - 1] * mat_radii_[i - 1];
    }
    // Get 1/8th the volume of the full ring, for each angular segment.
    vol /= 8.;

    // NNE
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
        pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
      fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // NEE
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
        pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
      fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // SEE
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
        pin_type_ == PinCellType::YN || pin_type_ == PinCellType::IV) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
      fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // SSE
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
        pin_type_ == PinCellType::YN || pin_type_ == PinCellType::IV) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
      fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // SSW
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
        pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
      fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // SWW
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
        pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
      fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // NWW
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
        pin_type_ == PinCellType::YP || pin_type_ == PinCellType::II) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
      fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }

    // NNW
    if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
        pin_type_ == PinCellType::YP || pin_type_ == PinCellType::II) {
      fsrs_.emplace_back();
      fsrs_.back().volume() = vol;
      fsrs_.back().xs() = mats_[i];
      fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
      fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
      fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
      if (i > 0) {
        fsrs_.back().tokens().push_back(
            {radii_[i - 1], Surface::Side::Positive});
      }
    }
  }

  // Must now make the 8 outer regions. We start by calculating their volumes.
  // If dx != dy, then the bits can have different volumes. Here is a diagram
  // for the case of dx > dy, on just the upper right corner.
  // y     dd
  // ^   _______
  // |   |  /: |
  // | d | / : |
  // |   |/__:_|
  // ------------> x
  const double vol_ang = (PI * mat_radii_.back() * mat_radii_.back()) / 8.;
  double px = dx;
  if (pin_type_ == PinCellType::XN || pin_type_ == PinCellType::XP ||
      pin_type_ == PinCellType::I || pin_type_ == PinCellType::II ||
      pin_type_ == PinCellType::III || pin_type_ == PinCellType::IV) {
    px = 2. * dx;
  }
  double py = dy;
  if (pin_type_ == PinCellType::YN || pin_type_ == PinCellType::YP ||
      pin_type_ == PinCellType::I || pin_type_ == PinCellType::II ||
      pin_type_ == PinCellType::III || pin_type_ == PinCellType::IV) {
    py = 2. * dy;
  }
  const double d = 0.5 * std::min(px, py);
  const double dd = 0.5 * std::max(px, py);
  const double vol_min = (0.5 * d * d) - vol_ang;
  const double vol_max = vol_min + (d * (dd - d));

  // NNE
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
      pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_min : vol_max;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});
  }

  // NEE
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
      pin_type_ == PinCellType::YP || pin_type_ == PinCellType::I) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_max : vol_min;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});
  }

  // SEE
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
      pin_type_ == PinCellType::YN || pin_type_ == PinCellType::IV) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_max : vol_min;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});
  }

  // SSE
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XP ||
      pin_type_ == PinCellType::YN || pin_type_ == PinCellType::IV) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_min : vol_max;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});
  }

  // SSW
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
      pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_min : vol_max;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});
  }

  // SWW
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
      pin_type_ == PinCellType::YN || pin_type_ == PinCellType::III) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_max : vol_min;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});
  }

  // NWW
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
      pin_type_ == PinCellType::YP || pin_type_ == PinCellType::II) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_max : vol_min;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});
  }

  // NNW
  if (pin_type_ == PinCellType::Full || pin_type_ == PinCellType::XN ||
      pin_type_ == PinCellType::YP || pin_type_ == PinCellType::II) {
    fsrs_.emplace_back();
    fsrs_.back().volume() = px > py ? vol_min : vol_max;
    fsrs_.back().xs() = mats_.back();
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});
  }
}

}  // namespace scarabee
