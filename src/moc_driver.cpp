#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>

#include <cmath>

MOCDriver::MOCDriver(std::shared_ptr<Cartesian2D> geometry, double d,
                     std::uint32_t na)
    : angle_info_(), tracks_(), geometry_(geometry), n_angles_(na) {
  // Determine the angles and spacings for the tracks
  double delta_phi = 2. * PI / static_cast<double>(n_angles_);

  // If we have reflective boundaries (which is currently the only option),
  // then we must divide n_angles_ by 2, as instead of needing all [0,2pi],
  // we only need [0, pi], and we can use tracks in the opposite direction.
  std::uint32_t n_track_angles_ = n_angles_ / 2;

  angle_info_.resize(n_track_angles_, {0., 0., 0., 0, 0});
  const double Dx = geometry_->x_max() - geometry_->x_min();
  const double Dy = geometry_->y_max() - geometry_->y_min();
  for (std::uint32_t i = 0; i < n_track_angles_; i++) {
    // Get the initial guess for phi_i
    double phi_i = delta_phi * (static_cast<double>(i) + 0.5);

    // Compute the number of x and y
    double nx = std::floor((Dy / d) * std::abs(std::cos(phi_i))) + 1.;
    double ny = std::floor((Dx / d) * std::abs(std::sin(phi_i))) + 1.;

    // Calculate information for a given angle, except the weight.
    // Weight is calculated once all angles are known.
    angle_info_[i].phi = std::atan((Dy * nx) / (Dx * ny));
    angle_info_[i].d = (Dx / nx) * std::sin(angle_info_[i].phi);
    angle_info_[i].nx = static_cast<std::uint32_t>(nx);
    angle_info_[i].ny = static_cast<std::uint32_t>(ny);
  }

  // Go through and calculate the angle weights
  for (std::uint32_t i = 0; i < n_track_angles_; i++) {
    if (i == 0) {
      const double phi_0 = angle_info_[0].phi;
      const double phi_1 = angle_info_[1].phi;
      angle_info_[i].wgt = (1. / (2. * PI)) * (0.5 * (phi_1 - phi_0) + phi_0);
    } else if (i == (n_angles_ - 1)) {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_i = angle_info_[i].phi;
      angle_info_[i].wgt =
          (1. / (2. * PI)) * (2. * PI - phi_i + 0.5 * (phi_i - phi_im1));
    } else {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_ip1 = angle_info_[i + 1].phi;
      angle_info_[i].wgt = (1. / (4. * PI)) * 0.5 * (phi_ip1 - phi_im1);
    }
  }

  // Now that we have all of the angle information, we can create tracks.
  tracks_.resize(n_track_angles_);
  for (std::uint32_t i = 0; i < n_track_angles_; i++) {
    const auto& ai = angle_info_[i];

    // Allocate space for tracks associated with this angle
    tracks_[i].reserve(ai.nx + ai.ny);

    // Get direction for the track
    Direction u(ai.phi);

    // spacing between starts in x
    const double dx = Dx / static_cast<double>(ai.nx);
    // spacing between starts in y
    const double dy = Dy / static_cast<double>(ai.ny);

    // Depending on angle, we either start on the -x bound, or the +y bound
    if (ai.phi < 0.5 * PI) {
      // Start on -x boundary in upper right corner and move down
      double x = geometry_->x_min();
      double y = dy * (static_cast<double>(ai.ny) + 0.5) + geometry_->y_min();
      for (std::uint32_t t = 0; t < (ai.nx + ai.ny); t++) {
        if (t == ai.ny) {
          // Next, we move across the -y boundary
          y = geometry_->y_min();
          x = 0.5 * dx + geometry_->x_min();
        }

        Vector r_start(x, y);
        Vector r_end = r_start;
        std::vector<Segment> segments;
        geometry_->trace_segments(r_end, u, segments);
        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt, segments);

        if (t < ai.ny) {
          y -= dy;
        } else {
          x += dx;
        }
      }
    } else {
      // Start on -y boundary in lower right corner and move across
      double x = geometry_->x_min();
      double y = 0.5 * dy + geometry_->y_min();
      for (std::uint32_t t = 0; t < (ai.nx + ai.ny); t++) {
        if (t == ai.nx) {
          // Next, we move across the -y boundary
          x = geometry_->x_max();
          y = 0.5 * dy + geometry_->y_min();
        }

        Vector r_start(x, y);
        Vector r_end = r_start;
        std::vector<Segment> segments;
        geometry_->trace_segments(r_end, u, segments);
        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt, segments);

        if (t < ai.nx) {
          x += dx;
        } else {
          y += dy;
        }
      }
    }
  }
}