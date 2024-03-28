#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xmath.hpp>

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
      flux_(),
      src_(),
      ngroups_(0),
      n_pol_angles_(polar_quad_.abscissae().size()),
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

  // Allocate arrays and assign indices
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    fsrs_[i]->set_indx(i);
  }

  flux_.resize({fsrs_.size(), ngroups_});
  src_.resize({fsrs_.size(), ngroups_});
}

void MOCDriver::set_flux_tolerance(double ftol) {
  if (ftol <= 0.) {
    throw ScarabeeException(
        "Tolerance for flux must be in the interval (0., 0.1).");
  }

  if (ftol >= 0.1) {
    throw ScarabeeException(
        "Tolerance for flux must be in the interval (0., 0.1).");
  }

  flux_tol_ = ftol;
}

void MOCDriver::set_keff_tolerance(double ktol) {
  if (ktol <= 0.) {
    throw ScarabeeException(
        "Tolerance for keff must be in the interval (0., 0.1).");
  }

  if (ktol >= 0.1) {
    throw ScarabeeException(
        "Tolerance for keff must be in the interval (0., 0.1).");
  }

  keff_tol_ = ktol;
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

  generate_azimuthal_quadrature(n_angles, d);
  generate_tracks();
  set_track_ends_bcs();
  allocate_track_fluxes();
  calculate_segment_exps();
}

void MOCDriver::solve_keff() {
  // Create a new array to hold the source according to the current flux, and
  // another for the next generation flux.
  xt::xtensor<double, 2> scat_source(flux_.shape());
  xt::xtensor<double, 2> fiss_source(flux_.shape());
  xt::xtensor<double, 2> next_flux(flux_.shape());
  flux_.fill(1.);
  next_flux.fill(1.);
  keff_ = calc_keff(flux_);
  double old_keff = 100.;
  std::size_t outer_iter = 0;

  // Initialize angular flux
  fill_fission_source(fiss_source, flux_);
  fill_scatter_source(scat_source, next_flux);
  src_ = fiss_source + scat_source;
  const double init_ang_flx = std::sqrt(xt::sum(src_*src_)());
  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      track.entry_flux().fill(init_ang_flx);
      track.exit_flux().fill(init_ang_flx);
    }
  }

  // Outer Generations
  while (std::abs(old_keff - keff_) > keff_tol_) {
    outer_iter++;

    // At the beginning of a generation, we calculate the fission source
    fill_fission_source(fiss_source, flux_);

    // Copy flux into next_flux, so that we can continually use next_flux
    // in the inner iterations for the scattering source.
    next_flux = flux_;

    double max_flux_diff = 100.;
    std::size_t inner_iter = 0;
    // Inner Iterations
    while (max_flux_diff > flux_tol_) {
      inner_iter++;

      // At the beginning of an inner iteration, we calculate the fission source
      fill_scatter_source(scat_source, next_flux);
      
      // Calculate the new source
      src_ = fiss_source + scat_source;

      next_flux.fill(0.);
      sweep(next_flux); 

      // Calculate the max difference in the flux
      max_flux_diff = 0.;
      for (std::size_t i = 0; i < fsrs_.size(); i++) {
        for (std::uint32_t g = 0; g < ngroups_; g++) {
          const double rel_diff = std::abs(next_flux(i, g) - flux_(i, g)) / flux_(i, g);
          if (rel_diff > max_flux_diff) max_flux_diff = rel_diff;
        }
      }

      // Copy next_flux into flux for calculating next relative difference
      flux_ = next_flux;

      std::cout << "     Finished inner iteration " << inner_iter << ", max_flux_diff = " << max_flux_diff << "\n";
      
      if (inner_iter > 100)
        break;
    } // End of Inner Iterations

    // Calculate keff
    old_keff = keff_;
    keff_ = calc_keff(next_flux);

    std::cout << " >> Iter " << outer_iter << " keff: " << keff_ << ", old_keff: " << old_keff << "\n";

    // Assign next_flux to be the flux
    std::swap(next_flux, flux_);
  }  // End of Outer Generations
}

void MOCDriver::sweep(xt::xtensor<double, 2>& sflux) {
  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      auto angflux = track.entry_flux();
      const double tw = track.weight();

      for (auto& seg : track) {
        const std::size_t i = seg.fsr_indx();
        const double seg_const = (4. * PI * seg.length() / seg.volume());

        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            const double delta_flx = (angflux(g, p) - src_(i, g) / seg.xs().Et(g)) * (1. - seg.exp()(g, p));

            if (delta_flx < 0.) {
              throw ScarabeeException("Neg delta_flux");
            } else if (seg_const < 0.) {
              std::cout.precision(16);
              std::cout << "Bad seg " << seg.length() << ", " << seg.volume() << "\n";
              throw ScarabeeException("Neg seg_const");
            } else if (tw < 0.) {
              throw ScarabeeException("Neg tw");
            } else if (polar_quad_.weights()[p] < 0.) {
              throw ScarabeeException("Neg pol weight");
            } else if (polar_quad_.abscissae()[p] < 0.) {
              throw ScarabeeException("Neg pol abs");
            }

            sflux(i, g) += seg_const * tw * polar_quad_.weights()[p] * polar_quad_.abscissae()[p] * delta_flx;
            angflux(g, p) -= delta_flx;
          }  // For all polar angles
        }    // For all groups

      }      // For all segments along forward direction of track

      // Set incoming flux for next track
      if (track.exit_bc() == BoundaryCondition::Reflective) {
        track.exit_track_flux() = angflux;
      } else {
        // Vacuum
        track.exit_track_flux().fill(0.);
      }

      // Do same track but in reverse
      angflux = track.exit_flux();
      for (auto seg_it = track.rbegin(); seg_it != track.rend(); seg_it++) {
        auto& seg = *seg_it;
        const std::size_t i = seg.fsr_indx();
        const double seg_const = (4. * PI * seg.length() / seg.volume());
        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            const double delta_flx = (angflux(g, p) - src_(i, g) / seg.xs().Et(g)) * (1. - seg.exp()(g, p));

            if (delta_flx < 0.) {
              throw ScarabeeException("Neg delta_flux");
            } else if (seg_const < 0.) {
              throw ScarabeeException("Neg seg_const");
            } else if (tw < 0.) {
              throw ScarabeeException("Neg tw");
            } else if (polar_quad_.weights()[p] < 0.) {
              throw ScarabeeException("Neg pol weight");
            } else if (polar_quad_.abscissae()[p] < 0.) {
              throw ScarabeeException("Neg pol abs");
            }
            sflux(i, g) += seg_const * tw * polar_quad_.weights()[p] * polar_quad_.abscissae()[p] * delta_flx;
            angflux(g, p) -= delta_flx;
          }  // For all polar angles
        }    // For all groups
      }      // For all segments along forward direction of track

      // Set incoming flux for next track
      if (track.entry_bc() == BoundaryCondition::Reflective) {
        track.entry_track_flux() = angflux;
      } else {
        // Vacuum
        track.entry_track_flux().fill(0.);
      }

    }  // For all tracks
  }    // For all azimuthal angles
}

double MOCDriver::calc_keff(const xt::xtensor<double, 2>& flux) const {
  double keff = 0.;
  
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const double Vr = fsrs_[i]->volume();
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      const double nu = mat.nu(g);
      const double Ef = mat.Ef(g);
      const double flx = flux(i, g);

      keff += Vr * nu * Ef * flx;
    }
  }

  return keff;
}

double MOCDriver::Qscat(std::uint32_t g, std::size_t i, const xt::xtensor<double, 2>& flux) const {
  double Qout = 0.;
  const auto& mat = *fsrs_[i]->xs();

  for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
    const double flux_gg_i = flux(gg, i);

    // Scattering into group g, excluding g -> g
    if (gg != g) Qout += mat.Es(gg, g) * flux_gg_i;
  }

  return Qout;
}

double MOCDriver::Qfiss(std::uint32_t g, std::size_t i, const xt::xtensor<double, 2>& flux) const {
  double Qout = 0.;
  const double inv_k = 1. / keff_;
  const auto& mat = *fsrs_[i]->xs();
  const double chi_g = mat.chi(g);

  for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
    const double Ef_gg = mat.Ef(gg);
    const double flux_gg_i = flux(gg, i);

    // Prompt Fission
    Qout += inv_k * chi_g * mat.nu(gg) * Ef_gg * flux_gg_i;
  }

  return Qout;
}

void MOCDriver::fill_scatter_source(xt::xtensor<double, 2>& scat_src, const xt::xtensor<double, 2>& flux) const {
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      scat_src(i, g) = Qscat(g, i, flux);
    }
  }
}

void MOCDriver::fill_fission_source(xt::xtensor<double, 2>& fiss_src, const xt::xtensor<double, 2>& flux) const {
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      fiss_src(i, g) = Qfiss(g, i, flux);
    }
  }
}

void MOCDriver::generate_azimuthal_quadrature(std::uint32_t n_angles,
                                              double d) {
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
      angle_info_[i].wgt = (1. / (2. * PI)) * (0.5 * (phi_1 - phi_0) + phi_0);
    } else if (i == (n_track_angles_ - 1)) {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_i = angle_info_[i].phi;
      angle_info_[i].wgt =
          (1. / (2. * PI)) * (PI - phi_i + 0.5 * (phi_i - phi_im1));
    } else {
      const double phi_im1 = angle_info_[i - 1].phi;
      const double phi_ip1 = angle_info_[i + 1].phi;
      angle_info_[i].wgt = (1. / (4. * PI)) * (phi_ip1 - phi_im1);
    }
  }
}

void MOCDriver::generate_tracks() {
  std::uint32_t n_track_angles_ =
      static_cast<std::uint32_t>(angle_info_.size());
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

    // Depending on angle, we either start on the -x bound, or the +x bound
    if (ai.phi < 0.5 * PI) {
      // Start on -x boundary in upper left corner and move down
      double x = geometry_->x_min();
      double y = dy * (static_cast<double>(ai.ny - 1) + 0.5) + geometry_->y_min();
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
        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt * ai.d,
                                segments);

        if (t < ai.ny) {
          y -= dy;
        } else {
          x += dx;
        }
      }
    } else {
      // Start on -y boundary in lower left corner and move across and up
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
        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt * ai.d,
                                segments);

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
    auto& comp_tracks = *(tracks_.rbegin() + a);

    // Go through intersections on top side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(i).set_exit_track_flux(&comp_tracks.at(ai.ny + i).exit_flux());
      comp_tracks.at(ai.ny + i).set_exit_track_flux(&tracks.at(i).exit_flux());

      tracks.at(i).exit_bc() = this->y_max_bc_;
      comp_tracks.at(ai.nx - 1 - i).exit_bc() = this->y_max_bc_;
    }

    // Go through intersections on bottom side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(ai.ny + i).set_entry_track_flux(&comp_tracks.at(i).entry_flux());
      comp_tracks.at(i).set_entry_track_flux(&tracks.at(ai.ny + i).entry_flux());

      tracks.at(ai.ny + i).entry_bc() = this->y_min_bc_;
      comp_tracks.at(nt - 1 - i).entry_bc() = this->y_min_bc_;
    }

    // Go down left/right sides
    for (std::uint32_t i = 0; i < ai.ny; i++) {
      // Left
      tracks.at(i).set_entry_track_flux(&comp_tracks.at(ai.ny - 1 - i).exit_flux());
      comp_tracks.at(ai.ny - 1 - i).set_exit_track_flux(&tracks.at(i).entry_flux());

      tracks.at(i).entry_bc() = this->x_min_bc_;
      comp_tracks.at(ai.nx + i).exit_bc() = this->x_min_bc_;

      // Right
      tracks.at(ai.nx + i).set_exit_track_flux(&comp_tracks.at(nt - 1 - i).entry_flux());
      comp_tracks.at(nt - 1 - i).set_entry_track_flux(&tracks.at(ai.nx + i).exit_flux());

      tracks.at(ai.nx + i).exit_bc() = this->x_max_bc_;
      comp_tracks.at(i).entry_bc() = this->x_max_bc_;
    }
  }
}

void MOCDriver::allocate_track_fluxes() {
  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      track.entry_flux().resize({ngroups_, n_pol_angles_});
      track.exit_flux().resize({ngroups_, n_pol_angles_});
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
        auto& exp = seg.exp();
        exp.resize({ngroups_, n_pol_angles});

        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles; p++) {
            exp(g, p) = std::exp(l * Et(g) * pd(p));
          }
        }

      }  // For all Segments
    }    // For all Tracks
  }      // For all angles
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
