#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

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
    auto mssg = "MOCDriver provided with nullptr geometry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (geometry_->tiles_valid() == false) {
    auto mssg = "Cannot run MOC on a Cartesian2D geometry with invalid tiles.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Get all FSRs
  auto nfsr = geometry_->num_fsrs();
  fsrs_.reserve(nfsr);
  geometry_->append_fsrs(fsrs_);

  if (fsrs_.empty()) {
    auto mssg = "No flat source regions found.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
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
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ftol >= 0.1) {
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  flux_tol_ = ftol;
}

void MOCDriver::set_keff_tolerance(double ktol) {
  if (ktol <= 0.) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ktol >= 0.1) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  keff_tol_ = ktol;
}

void MOCDriver::draw_tracks(std::uint32_t n_angles, double d) {
  if (n_angles % 2 != 0) {
    // If the number of angles is odd, an angle will be lost
    auto mssg = "MOCDriver must have an even number of angles.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (n_angles < 4) {
    auto mssg = "MOCDriver must have at least 4 angles.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (d <= 0.) {
    auto mssg = "MOCDriver track spacing must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Clear any previous data
  angle_info_.clear();
  tracks_.clear();

  generate_azimuthal_quadrature(n_angles, d);
  generate_tracks();
  segment_renormalization();
  calculate_segment_exps();
  set_track_ends_bcs();
  allocate_track_fluxes();
}

void MOCDriver::solve_keff() {
  spdlog::info("Solving for keff.");
  spdlog::info("keff tolerance: {:.5E}", keff_tol_);
  spdlog::info("Flux tolerance: {:.5E}", flux_tol_);

  flux_.resize({fsrs_.size(), ngroups_});
  src_.resize({fsrs_.size(), ngroups_});
  src_.fill(0.);
  auto fiss_src = src_;
  auto self_scat_src = src_;
  auto prev_self_scat_src = src_;
  auto extern_scat_src = src_;
  auto extern_src = src_;

  // Initialize flux and keff
  flux_.fill(1.);
  keff_ = 1.;
  auto next_flux = flux_;
  double prev_keff = keff_;

  spdlog::info("Initial keff {:.5f}", keff_);

  // Initialize angular flux
  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      track.entry_flux().fill(1./(4.*PI));
      track.exit_flux().fill(1./(4.*PI));
    }
  }

  double rel_diff_keff = 100.;
  double max_flx_diff = 100;
  std::size_t outer_iter = 0;
  while (rel_diff_keff > keff_tol_ || max_flx_diff > flux_tol_) {
    outer_iter++;

    fill_fission_source(fiss_src, next_flux);
    fill_external_scatter_source(extern_scat_src, next_flux);
    fill_self_scatter_source(self_scat_src, next_flux);
    extern_src = fiss_src + extern_scat_src;

    double max_src_diff = 100.;
    std::size_t inner_iter = 0;
    while (max_src_diff > flux_tol_) {
      inner_iter++;

      src_ = extern_src + self_scat_src;

      sweep(next_flux);

      // Calculate new self source for next inner iteration
      prev_self_scat_src = self_scat_src;
      fill_self_scatter_source(self_scat_src, next_flux);

      // Get difference 
      max_src_diff = xt::amax(xt::abs(self_scat_src - prev_self_scat_src) / self_scat_src)();

      spdlog::debug("Inner iteration {} source diff {:.5E}", inner_iter, max_src_diff);
    }

    prev_keff = keff_;
    keff_ = calc_keff(next_flux, flux_);
    rel_diff_keff = std::abs(keff_ - prev_keff) / keff_;

    // Get difference 
    max_flx_diff = xt::amax(xt::abs(next_flux - flux_) / next_flux)();

    flux_ = next_flux;
    spdlog::info("Iteration {} keff {:.5f} flux difference {:.5E}", outer_iter, keff_, max_flx_diff);
  }
}

void MOCDriver::sweep(xt::xtensor<double, 2>& sflux) {
  for (auto& tracks : tracks_) {
#ifdef SCARABEE_USE_OMP
#pragma omp parallel for
#endif
    for (std::size_t t = 0; t < tracks.size(); t++) {
      auto& track = tracks[t];
      auto angflux = track.entry_flux();
      const double tw = track.weight();  // Azimuthal weight * track width

      // Follow track in forward direction
      for (auto& seg : track) {
        const std::size_t i = seg.fsr_indx();
        const double seg_const = (4. * PI * seg.length() / seg.volume());
        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            const double delta_flx = (angflux(g, p) - src_(i, g) / seg.xs().Et(g)) * (1. - seg.exp()(g, p));
#ifdef SCARABEE_USE_OMP
#pragma omp atomic
#endif
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

      // Follow track in backwards direction
      angflux = track.exit_flux();
      for (auto seg_it = track.rbegin(); seg_it != track.rend(); seg_it++) {
        auto& seg = *seg_it;
        const std::size_t i = seg.fsr_indx();
        const double seg_const = (4. * PI * seg.length() / seg.volume());
        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            const double delta_flx = (angflux(g, p) - src_(i, g) / seg.xs().Et(g)) * (1. - seg.exp()(g, p));
#ifdef SCARABEE_USE_OMP
#pragma omp atomic
#endif
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

/*
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
*/

double MOCDriver::calc_keff(const xt::xtensor<double, 2>& flux, const xt::xtensor<double, 2>& old_flux) const {
  double num = 0.;
  double denom = 0.;

  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const double Vr = fsrs_[i]->volume();
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      const double nu = mat.nu(g);
      const double Ef = mat.Ef(g);
      const double VvEf = Vr * nu * Ef;
      const double flx = flux(i, g);
      const double oflx = old_flux(i, g);

      if (flx == oflx) {
        spdlog::error("New and old flux is equal");
      }

      num += VvEf * flx;
      denom += VvEf * oflx;
    }
  }

  return keff_ * num / denom;
}

/*
double MOCDriver::calc_abs(const xt::xtensor<double, 2>& flux) const {
  double abs = 0.;

  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const double Vr = fsrs_[i]->volume();
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      const double Ea = mat.Ea(g);
      const double flx = flux(i, g);
      abs += Vr * Ea * flx;
    }
  }

  return abs;
}
*/

/*
void MOCDriver::fill_scatter_source(xt::xtensor<double, 2>& scat_src,
                                    const xt::xtensor<double, 2>& flux) const {
  const double isotropic = 1. / (4. * PI);
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      double Qout = 0.;

      for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
        const double flux_gg_i = flux(i, gg);
        const double Es_gg_to_g = mat.Es(gg, g);
        Qout += Es_gg_to_g * flux_gg_i;
      }

      scat_src(i, g) = isotropic * Qout;
    }
  }
}
*/

void MOCDriver::fill_external_scatter_source(xt::xtensor<double, 2>& scat_src, const xt::xtensor<double, 2>& flux) const {
  const double isotropic = 1. / (4. * PI);
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      double Qout = 0.;

      for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
        if (gg == g) continue; // Don't include self scattering
        const double flux_gg_i = flux(i, gg);
        const double Es_gg_to_g = mat.Es(gg, g);
        Qout += Es_gg_to_g * flux_gg_i;
      }

      scat_src(i, g) = isotropic * Qout;
    }
  }
}

void MOCDriver::fill_self_scatter_source(xt::xtensor<double, 2>& scat_src, const xt::xtensor<double, 2>& flux) const {
  const double isotropic = 1. / (4. * PI);
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      scat_src(i, g) = isotropic * mat.Es(g, g) * flux(i, g);
    }
  }
}

void MOCDriver::fill_fission_source(xt::xtensor<double, 2>& fiss_src,
                                    const xt::xtensor<double, 2>& flux) const {
  const double inv_k = 1. / keff_;
  const double isotropic = 1. / (4. * PI);
  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const auto& mat = *fsrs_[i]->xs();

    double fiss_rate_i = 0;
    for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
      const double nu_gg = mat.nu(gg);
      const double Ef_gg = mat.Ef(gg);
      const double flux_gg_i = flux(i, gg);
      fiss_rate_i += nu_gg * Ef_gg * flux_gg_i;
    }

    for (std::uint32_t g = 0; g < ngroups_; g++) {
      const double chi_g = mat.chi(g);
      fiss_src(i, g) = isotropic * inv_k * chi_g * fiss_rate_i;
    }
  }
}

void MOCDriver::generate_azimuthal_quadrature(std::uint32_t n_angles,
                                              double d) {
  spdlog::info("Creating quadrature");

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

    const auto& ai = angle_info_[i];
    spdlog::debug("Angle {:.5E} pi: weight {:.5E}, width {:.5E}, nx {}, ny {}", ai.phi/PI, ai.wgt, ai.d, ai.nx, ai.ny);
  }
}

void MOCDriver::generate_tracks() {
  spdlog::info("Tracing tracks");

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
  spdlog::info("Determining track connections");
  
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

      if (tracks.at(i).exit_pos() != comp_tracks.at(ai.ny+i).exit_pos()) {
        std::stringstream mssg;
        mssg << "Disagreement in track end alignments: " << tracks.at(i).exit_pos() << " and " << comp_tracks.at(ai.ny+i).exit_pos() << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }

      tracks.at(i).exit_bc() = this->y_max_bc_;
      comp_tracks.at(ai.nx - 1 - i).exit_bc() = this->y_max_bc_;
    }

    // Go through intersections on bottom side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(ai.ny + i).set_entry_track_flux(
          &comp_tracks.at(i).entry_flux());
      comp_tracks.at(i).set_entry_track_flux(
          &tracks.at(ai.ny + i).entry_flux());

      if (tracks.at(ai.ny+i).entry_pos() != comp_tracks.at(i).entry_pos()) {
        std::stringstream mssg;
        mssg << "Disagreement in track end alignments: " << tracks.at(ai.ny+i).entry_pos() << " and " << comp_tracks.at(i).entry_pos() << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }

      tracks.at(ai.ny + i).entry_bc() = this->y_min_bc_;
      comp_tracks.at(nt - 1 - i).entry_bc() = this->y_min_bc_;
    }

    // Go down left/right sides
    for (std::uint32_t i = 0; i < ai.ny; i++) {
      // Left
      tracks.at(i).set_entry_track_flux(
          &comp_tracks.at(ai.ny - 1 - i).exit_flux());
      comp_tracks.at(ai.ny - 1 - i)
          .set_exit_track_flux(&tracks.at(i).entry_flux());
      
      if (tracks.at(i).entry_pos() != comp_tracks.at(ai.ny-1-i).exit_pos()) {
        std::stringstream mssg;
        mssg << "Disagreement in track end alignments: " << tracks.at(i).entry_pos() << " and " << comp_tracks.at(ai.ny-1-i).exit_pos() << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }

      tracks.at(i).entry_bc() = this->x_min_bc_;
      comp_tracks.at(ai.nx + i).exit_bc() = this->x_min_bc_;

      // Right
      tracks.at(ai.nx + i).set_exit_track_flux(
          &comp_tracks.at(nt - 1 - i).entry_flux());
      comp_tracks.at(nt - 1 - i)
          .set_entry_track_flux(&tracks.at(ai.nx + i).exit_flux());

      if (tracks.at(ai.nx+i).exit_pos() != comp_tracks.at(nt-1-i).entry_pos()) {
        std::stringstream mssg;
        mssg << "Disagreement in track end alignments: " << tracks.at(ai.nx+i).exit_pos() << " and " << comp_tracks.at(nt-1-i).entry_pos() << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }

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

void MOCDriver::segment_renormalization() {
  spdlog::info("Renormalizing segment lengths");

  // We now bias the traced segment lengths, so that we better predict the
  // volume of our flat source regions. We do this for each angle, but it
  // could be done for all angles together. A great explanation of this is
  // found in the MPACT theory manual ORNL/SPR-2021/2330 end of 5.4.

  // This holds the approximations for the FSR areas
  std::vector<double> approx_vols(fsrs_.size(), 0.);

  // Go through all angles
  for (std::size_t a = 0; a < angle_info_.size(); a++) {
    const double d = angle_info_[a].d; // Track width
    auto& tracks = tracks_[a]; // Vector of tracks

    // Zero approximate volumes
    for (auto& av : approx_vols) av = 0.;

    // Iterate through tracks and segments, adding contributions to volumes
    for (auto& track : tracks) {
      for (auto& seg : track) {
        const std::size_t i = seg.fsr_indx();
        approx_vols[i] += seg.length() * d;
      }
    }

    // Now we apply the corrections to the segment lengths 
    for (auto& track : tracks) {
      for (auto& seg : track) {
        const std::size_t i = seg.fsr_indx();
        seg.set_length(seg.length() * seg.volume() / approx_vols[i]);
      }
    }
  }
}

void MOCDriver::calculate_segment_exps() {
  spdlog::info("Calculating segment exponentials");

  xt::xtensor<double, 1> pd;  // Temp array to hold polar angle distance factors
  pd.resize({n_pol_angles_});
  for (std::size_t pi = 0; pi < n_pol_angles_; pi++) {
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
        exp.resize({ngroups_, n_pol_angles_});

        for (std::size_t g = 0; g < ngroups_; g++) {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            exp(g, p) = std::exp(l * Et(g) * pd(p));
          }  // For all polar angles
        }    // For all groups
      }      // For all Segments
    }        // For all Tracks
  }          // For all angles
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
