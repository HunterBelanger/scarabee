#include <moc/pin_cell.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>

#include <algorithm>

PinCell::PinCell(const std::vector<double>& mat_rads,
                 const std::vector<std::shared_ptr<TransportXS>>& mats,
                 std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
                 std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax)
    : Cell(xmin, xmax, ymin, ymax),
      mat_radii_(mat_rads),
      mats_(mats),
      radii_(),
      xm_(),
      pd_(),
      ym_(),
      nd_(),
      x_orig_(),
      y_orig_(),
      x0_(),
      y0_() {
  x_orig_ = 0.5 * (x_min_->x0() + x_max_->x0());
  y_orig_ = 0.5 * (y_min_->y0() + y_max_->y0());
  const double dx = (x_max_->x0() - x_min_->x0());
  const double dy = (y_max_->y0() - y_min_->y0());

  // For now, we center rings at the cell origin
  x0_ = x_orig_;
  y0_ = y_orig_;

  // Create the 4 surfaces which give us our 8 angular sections
  xm_ = std::make_shared<Surface>();
  xm_->type() = Surface::Type::XPlane;
  xm_->x0() = x0_;

  ym_ = std::make_shared<Surface>();
  ym_->type() = Surface::Type::YPlane;
  ym_->y0() = y0_;

  pd_ = std::make_shared<Surface>();
  pd_->type() = Surface::Type::Plane;
  pd_->A() = -1;
  pd_->B() = 1.;
  pd_->C() = y0_ - x0_;

  nd_ = std::make_shared<Surface>();
  nd_->type() = Surface::Type::Plane;
  nd_->A() = 1;
  nd_->B() = 1.;
  nd_->C() = y0_ + x0_;

  // Make sure the numbers for mats and rads is coherent
  if (mat_radii_.size() + 1 != mats_.size()) {
    throw ScarabeeException("Must have one more material than material radii.");
  }

  // Make sure radii all > 0
  for (auto r : mat_radii_) {
    if (r <= 0.) {
      throw ScarabeeException("All radii must be > 0.");
    }
  }

  // Make sure radii are sorted
  if (std::is_sorted(mat_radii_.begin(), mat_radii_.end()) == false) {
    throw ScarabeeException("All radii must be sorted.");
  }

  // Make sure we have at least 2 mats and one radii
  if (mat_radii_.size() == 0) {
    throw ScarabeeException("Must have at least one radiius.");
  }

  // Make sure largest radiius doesn't intersect the cell walls
  if (x0_ - mat_radii_.back() <= x_min_->x0() ||
      x0_ + mat_radii_.back() >= x_max_->x0() ||
      y0_ - mat_radii_.back() <= y_min_->y0() ||
      y0_ + mat_radii_.back() >= y_max_->y0()) {
    throw ScarabeeException(
        "Outer material radiius may not intersect cell wall.");
  }

  // First, make all anular regions
  for (std::size_t i = 0; i < mat_radii_.size(); i++) {
    double r = mat_radii_[i];
    radii_.push_back(std::make_shared<Surface>());
    radii_.back()->type() = Surface::Type::Cylinder;
    radii_.back()->x0() = x0_;
    radii_.back()->y0() = y0_;
    radii_.back()->r() = r;

    // Make 8 FSRs
    double vol = PI * r * r;
    if (i > 0) {
      vol -= PI * mat_radii_[i - 1] * mat_radii_[i - 1];
    }
    vol /= 8.;

    // NNE
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // NEE
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // SEE
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // SSE
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // SSW
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // SWW
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // NWW
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }

    // NNW
    fsrs_.emplace_back();
    fsrs_.back().volume() = vol;
    fsrs_.back().xs() = mats_[i];
    fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Negative});
    fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
    fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
    if (i > 0) {
      fsrs_.back().tokens().push_back({radii_[i - 1], Surface::Side::Positive});
    }
  }

  // Must now make the 8 outer regions
  double vol = ((dx * dy) - (PI * mat_radii_.back() * mat_radii_.back())) / 8.;

  // NNE
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});

  // NEE
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});

  // SEE
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({x_max_, Surface::Side::Negative});

  // SSE
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({xm_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});

  // SSW
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({pd_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({y_min_, Surface::Side::Positive});

  // SWW
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({ym_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({pd_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});

  // NWW
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({ym_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({nd_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({x_min_, Surface::Side::Positive});

  // NNW
  fsrs_.emplace_back();
  fsrs_.back().volume() = vol;
  fsrs_.back().xs() = mats.back();
  fsrs_.back().tokens().push_back({radii_.back(), Surface::Side::Positive});
  fsrs_.back().tokens().push_back({xm_, Surface::Side::Negative});
  fsrs_.back().tokens().push_back({nd_, Surface::Side::Positive});
  fsrs_.back().tokens().push_back({y_max_, Surface::Side::Negative});
}