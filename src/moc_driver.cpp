#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>

#include <cmath>

MOCDriver::MOCDriver(std::shared_ptr<Cartesian2D> geometry,
                     PolarQuadrature polar_quad, BoundaryCondition xmin,
                     BoundaryCondition xmax, BoundaryCondition ymin,
                     BoundaryCondition ymax)
    : angle_info_(),
      tracks_(),
      fsrs_(),
      geometry_(geometry),
      polar_quad_(polar_quad),
      ngroups_(0),
      x_min_bc_(xmin),
      x_max_bc_(xmax),
      y_min_bc_(ymin),
      y_max_bc_(ymax) {
  if (geometry_ == nullptr) {
    throw ScarabeeException("MOCDriver provided with nullptr geometry.");
  }

  if (geometry_->tiles_valid() == false) {
    throw ScarabeeException(
        "Cannot run MOC on a Cartesian2D geometry with invalid tiles.");
  }

  // Get all FSRs
  auto nfsr = geometry_->num_fsrs();
  fsrs_.reserve(nfsr);
  geometry_->append_fsrs(fsrs_);

  if (fsrs_.empty()) {
    throw ScarabeeException("No flat source regions found.");
  }

  ngroups_ = fsrs_.front()->xs()->ngroups();
}

void MOCDriver::draw_tracks(std::uint32_t n_angles, double d) {
  if (n_angles % 2 != 0) {
    // If the number of angles is odd, an angle will be lost
    throw ScarabeeException("MOCDriver must have an even number of angles.");
  }

  if (n_angles < 4) {
    throw ScarabeeException("MOCDriver must have at least 4 angles.");
  }

  if (d <= 0.) {
    throw ScarabeeException("MOCDriver track spacing must be > 0.");
  }

  // Clear any previous data
  angle_info_.clear();
  tracks_.clear();

  generate_azimuthal_quadrature(n_angles, d) ;
  generate_tracks();
  set_track_ends_bcs();
  allocate_track_fluxes();
  allocate_fsr_flux_source();
  calculate_segment_exps();
}

void MOCDriver::generate_azimuthal_quadrature(std::uint32_t n_angles, double d) {
  // Determine the angles and spacings for the tracks
  double delta_phi = 2. * PI / static_cast<double>(n_angles);

  // If we have reflective boundaries (which is currently the only option),
  // then we must divide n_angles by 2, as instead of needing all [0,2pi],
  // we only need [0, pi], and we can use tracks in the opposite direction.
  std::uint32_t n_track_angles_ = n_angles / 2;

  angle_info_.resize(n_track_angles_, {0., 0., 0., 0, 0});
  const double Dx = geometry_->x_max() - geometry_->x_min();
  const double Dy = geometry_->y_max() - geometry_->y_min();
  for (std::uint32_t i = 0; i < n_track_angles_; i++) {
    // Get the initial guess for phi_i
    double phi_i = delta_phi * (static_cast<double>(i) + 0.5);

    // Compute the number of x and y
    double nx = std::floor((Dy / d) * std::abs(std::sin(phi_i))) + 1.;
    double ny = std::floor((Dx / d) * std::abs(std::cos(phi_i))) + 1.;

    // Calculate information for a given angle, except the weight.
    // Weight is calculated once all angles are known.
    angle_info_[i].phi = std::atan((Dy * nx) / (Dx * ny));

    // Fix angles between [pi/2, pi], due to arctan result domain
    if (phi_i > 0.5 * PI) angle_info_[i].phi = PI - angle_info_[i].phi;

    angle_info_[i].d = (Dx / nx) * std::sin(angle_info_[i].phi);
    angle_info_[i].nx = static_cast<std::uint32_t>(nx);
    angle_info_[i].ny = static_cast<std::uint32_t>(ny);
  }

  // Go through and calculate the angle weights
  for (std::uint32_t i = 0; i < n_track_angles_; i++) {
    if (i == 0) {
      const double phi_0 = angle_info_[0].phi;
      const double phi_1 = angle_info_[1].phi;
      angle_info_[i].wgt = (2. / (2. * PI)) * (0.5 * (phi_1 - phi_0) + phi_0);
    } else if (i == (n_track_angles_ - 1)) {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_i = angle_info_[i].phi;
      angle_info_[i].wgt =
          (2. / (2. * PI)) * (PI - phi_i + 0.5 * (phi_i - phi_im1));
    } else {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_ip1 = angle_info_[i + 1].phi;
      angle_info_[i].wgt = (2. / (4. * PI)) * (phi_ip1 - phi_im1);
    }
  }
}

void MOCDriver::generate_tracks() {
  std::uint32_t n_track_angles_ = static_cast<std::uint32_t>(angle_info_.size());
  const double Dx = geometry_->x_max() - geometry_->x_min();
  const double Dy = geometry_->y_max() - geometry_->y_min();

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
      double y =
          dy * (static_cast<double>(ai.ny - 1) + 0.5) + geometry_->y_min();
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
      double x = 0.5 * dx + geometry_->x_min();
      double y = geometry_->y_min();
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

void MOCDriver::set_track_ends_bcs() {
  // We only go through the first half of the tracks, where phi < pi / 2.
  for (std::size_t a = 0; a < angle_info_.size() / 2; a++) {
    const auto& ai = angle_info_[a];
    const std::uint32_t nt = ai.nx + ai.ny;
    auto& tracks = tracks_[a];
    auto& comp_tracks = tracks_[angle_info_.size() - 1 - a];

    // Go through intersections on top side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(i).set_exit_track(&comp_tracks.at(ai.nx - 1 - i));
      comp_tracks.at(ai.nx - 1 - i).set_exit_track(&tracks.at(i));

      tracks.at(i).exit_bc() = this->y_max_bc_;
      comp_tracks.at(ai.nx - 1 - i).exit_bc() = this->y_max_bc_;
    }

    // Go through intersections on bottom side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(ai.ny + i).set_entry_track(&comp_tracks.at(nt - 1 - i));
      comp_tracks.at(nt - 1 - i).set_entry_track(&tracks.at(ai.ny + i));

      tracks.at(ai.ny + i).entry_bc() = this->y_min_bc_;
      comp_tracks.at(nt - 1 - i).entry_bc() = this->y_min_bc_;
    }

    // Go down left/right sides
    for (std::uint32_t i = 0; i < ai.ny; i++) {
      // Left
      tracks.at(i).set_entry_track(&comp_tracks.at(ai.nx + i));
      comp_tracks.at(ai.nx + i).set_exit_track(&tracks.at(i));

      tracks.at(i).entry_bc() = this->x_min_bc_;
      comp_tracks.at(ai.nx + i).exit_bc() = this->x_min_bc_;

      // Right
      tracks.at(ai.nx + i).set_exit_track(&comp_tracks.at(i));
      comp_tracks.at(i).set_entry_track(&tracks.at(ai.nx + i));

      tracks.at(ai.nx + i).exit_bc() = this->x_max_bc_;
      comp_tracks.at(i).entry_bc() = this->x_max_bc_;
    }
  }
}

void MOCDriver::allocate_track_fluxes() {
  // Get the number of groups, and number of polar angles
  std::size_t n_pol_angles = polar_quad_.abscissae().size();

  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      track.entry_flux().reshape({ngroups_, n_pol_angles});
      track.exit_flux().reshape({ngroups_, n_pol_angles});
    }
  }
}

void MOCDriver::calculate_segment_exps() {
  // Get the number of groups, and number of polar angles
  const std::size_t n_pol_angles = polar_quad_.abscissae().size();
  xt::xtensor<double, 1> pd;
  pd.resize({n_pol_angles});
  for (std::size_t pi = 0; pi < n_pol_angles; pi++) {
    pd(pi) = 1. / polar_quad_.abscissae()[pi];
  }

  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      for (auto& seg : track) {
        const auto& Et = seg.xs().Et();
        const double l = -seg.length();
        // We unfortunately can't use xtensor-blas, as we don't have BLAS or
        // LAPACK on Windows. We instead construct the matrix ourselves.
        seg.exp().resize({ngroups_, n_pol_angles});

        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles; p++) {
            seg.exp()(g,p) = std::exp(l * Et(g) * pd(p));
          }
        }

      } // For all Segments
    } // For all Tracks
  } // For all angles
}

void MOCDriver::allocate_fsr_flux_source() {
  for (auto fsr : fsrs_) {
    fsr->flux().reshape({ngroups_});
    fsr->source().reshape({ngroups_});
  }
}

FlatSourceRegion& MOCDriver::get_fsr(const Vector& r, const Direction& u) {
  try {
    return geometry_->get_fsr(r, u);
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not find FSR in Cartesian2D.");
    throw err;
  }
}

const FlatSourceRegion& MOCDriver::get_fsr(const Vector& r,
                                           const Direction& u) const {
  try {
    return geometry_->get_fsr(r, u);
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not find FSR in Cartesian2D.");
    throw err;
  }
}