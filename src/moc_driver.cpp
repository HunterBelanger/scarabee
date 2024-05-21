#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>
#include <utils/timer.hpp>
#include <utils/math.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

#include <cmath>

namespace scarabee {

MOCDriver::MOCDriver(std::shared_ptr<Cartesian2D> geometry,
                     BoundaryCondition xmin, BoundaryCondition xmax,
                     BoundaryCondition ymin, BoundaryCondition ymax)
    : angle_info_(),
      tracks_(),
      geometry_(geometry),
      polar_quad_(YamamotoTabuchi<6>()),
      flux_(),
      extern_src_(),
      ngroups_(0),
      nfsrs_(0),
      n_pol_angles_(polar_quad_.sin().size()),
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

  // Get total number of unique FSRs
  nfsrs_ = geometry_->num_fsrs();
  if (nfsrs_ == 0) {
    auto mssg = "No flat source regions found.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  fsrs_.reserve(nfsrs_);

  // Get all FSRs IDs, and pointers to them
  auto fsr_ids = geometry_->get_all_fsr_ids();
  std::map<std::size_t, const FlatSourceRegion*> fsr_ptrs;
  geometry_->fill_fsrs(fsr_ptrs);

  // We now create offsets for each FSR ID, so that we can get a linear index
  // in the global FSR array, based on the ID and the unique instance number.
  for (auto id : fsr_ids) fsr_offsets_[id] = 0;
  auto id_it_prev = fsr_offsets_.begin();
  for (auto id_it = ++fsr_offsets_.begin(); id_it != fsr_offsets_.end();
       id_it++) {
    const std::size_t id_prev = id_it_prev->first;
    const std::size_t ninst_prev = geometry_->get_num_fsr_instances(id_prev);
    fsr_offsets_[id_it->first] = fsr_offsets_[id_prev] + ninst_prev;
    id_it_prev++;

    // Save ninst_prev points
    for (std::size_t i = 0; i < ninst_prev; i++) {
      fsrs_.push_back(fsr_ptrs[id_prev]);
    }
  }
  // Save pointers for last FSR
  const std::size_t id_prev = id_it_prev->first;
  const std::size_t ninst_prev = geometry_->get_num_fsr_instances(id_prev);
  for (std::size_t i = 0; i < ninst_prev; i++) {
    fsrs_.push_back(fsr_ptrs[id_prev]);
  }

  ngroups_ = geometry_->ngroups();

  // Allocate arrays and assign indices
  flux_.resize({ngroups_, nfsrs_});
  flux_.fill(0.);
  extern_src_.resize({ngroups_, nfsrs_});
  extern_src_.fill(0.);
}

std::size_t MOCDriver::size() const { return fsrs_.size(); }

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

void MOCDriver::generate_tracks(std::uint32_t n_angles, double d,
                                PolarQuadrature polar_quad) {
  // Timer for method
  Timer draw_timer;
  draw_timer.start();

  // Since we are reallocating bits, set solved to false
  solved_ = false;

  polar_quad_ = polar_quad;
  n_pol_angles_ = polar_quad_.sin().size();

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
  set_track_ends_bcs();
  allocate_track_fluxes();

  draw_timer.stop();
  spdlog::info("Time spent dawing tracks: {:.5} s.", draw_timer.elapsed_time());
}

void MOCDriver::solve() {
  Timer sim_timer;
  sim_timer.start();

  // Make sure the geometry has been drawn
  if (angle_info_.empty()) {
    auto mssg = "Cannot solve MOC problem. Geometry has not been traced.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (mode_ == SimulationMode::Keff) {
    spdlog::info("Solving for keff.");
    spdlog::info("keff tolerance: {:.5E}", keff_tol_);
  } else if (mode_ == SimulationMode::FixedSource) {
    spdlog::info("Solving fixed source problem.");
  }
  spdlog::info("Flux tolerance: {:.5E}", flux_tol_);

  if (mode_ == SimulationMode::FixedSource) {
    // Make sure extern_src_ is not all zero. Otherwise, we have no source !
    bool all_zero_extern_src = true;
    for (const auto& es : extern_src_) {
      if (es < 0.) {
        auto mssg = "All external sources must be > 0.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }

      if (es > 0.) {
        all_zero_extern_src = false;
      }
    }

    if (all_zero_extern_src) {
      auto mssg = "Must have at least one non-zero external source.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  flux_.resize({ngroups_, nfsrs_});
  xt::xtensor<double, 2> src_;
  src_.resize({ngroups_, nfsrs_});
  src_.fill(0.);

  // Initialize flux and keff
  if (mode_ == SimulationMode::Keff) {
    flux_.fill(1.);
  } else {
    flux_.fill(0.);
  }
  keff_ = 1.;
  auto next_flux = flux_;
  double prev_keff = keff_;

  // Initialize angular flux
  for (auto& tracks : tracks_) {
    for (auto& track : tracks) {
      track.entry_flux().fill(1. / (4. * PI));
      track.exit_flux().fill(1. / (4. * PI));
    }
  }

  double rel_diff_keff = 100.;
  if (mode_ == SimulationMode::FixedSource) {
    rel_diff_keff = 0.;
  }
  double max_flx_diff = 100;
  std::size_t iteration = 0;
  Timer iteration_timer;
  while (rel_diff_keff > keff_tol_ || max_flx_diff > flux_tol_) {
    iteration_timer.reset();
    iteration_timer.start();
    iteration++;

    fill_source(src_, flux_);
    src_ += extern_src_;

    next_flux.fill(0.);
    sweep(next_flux, src_);

    if (mode_ == SimulationMode::Keff) {
      prev_keff = keff_;
      keff_ = calc_keff(next_flux, flux_);
      rel_diff_keff = std::abs(keff_ - prev_keff) / keff_;
    }

    // Get difference
    max_flx_diff = xt::amax(xt::abs(next_flux - flux_) / next_flux)();

    flux_ = next_flux;

    iteration_timer.stop();
    spdlog::info("-------------------------------------");
    if (mode_ == SimulationMode::Keff) {
      spdlog::info("Iteration {:>4d}          keff: {:.5f}", iteration, keff_);
      spdlog::info("     keff difference:     {:.5E}", rel_diff_keff);
    } else if (mode_ == SimulationMode::FixedSource) {
      spdlog::info("Iteration {:>4d}", iteration);
    }
    spdlog::info("     max flux difference: {:.5E}", max_flx_diff);
    spdlog::info("     Iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());
  }

  solved_ = true;

  sim_timer.stop();
  spdlog::info("");
  spdlog::info("Simulation Time: {:.5E} s", sim_timer.elapsed_time());
}

void MOCDriver::sweep(xt::xtensor<double, 2>& sflux,
                      const xt::xtensor<double, 2>& src) {
#pragma omp parallel for
  for (std::size_t g = 0; g < ngroups_; g++) {
    for (auto& tracks : tracks_) {
      for (std::size_t t = 0; t < tracks.size(); t++) {
        auto& track = tracks[t];
        auto angflux = xt::view(track.entry_flux(), g, xt::all());
        const double tw = track.weight();  // Azimuthal weight * track width

        // Follow track in forward direction
        for (auto& seg : track) {
          const std::size_t i = seg.fsr_indx();
          const double Et = seg.xs()->Etr(g);
          const double Q = src(g, i);
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            double exp_m1 = mexp(Et * seg.length() * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux(p) - (Q / Et)) * exp_m1;
            sflux(g, i) += 4. * PI * tw * polar_quad_.wgt()[p] *
                           polar_quad_.sin()[p] * delta_flx;
            angflux(p) -= delta_flx;
          }  // For all polar angles
        }    // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.exit_bc() == BoundaryCondition::Reflective) {
          xt::view(track.exit_track_flux(), g, xt::all()) = angflux;
        } else {
          // Vacuum
          xt::view(track.exit_track_flux(), g, xt::all()).fill(0.);
        }

        // Follow track in backwards direction
        angflux = xt::view(track.exit_flux(), g, xt::all());
        for (auto seg_it = track.rbegin(); seg_it != track.rend(); seg_it++) {
          auto& seg = *seg_it;
          const std::size_t i = seg.fsr_indx();
          const double Et = seg.xs()->Etr(g);
          const double Q = src(g, i);
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            double exp_m1 = mexp(Et * seg.length() * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux(p) - (Q / Et)) * exp_m1;
            sflux(g, i) += 4. * PI * tw * polar_quad_.wgt()[p] *
                           polar_quad_.sin()[p] * delta_flx;
            angflux(p) -= delta_flx;
          }  // For all polar angles
        }    // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.entry_bc() == BoundaryCondition::Reflective) {
          xt::view(track.entry_track_flux(), g, xt::all()) = angflux;
        } else {
          // Vacuum
          xt::view(track.entry_track_flux(), g, xt::all()).fill(0.);
        }
      }  // For all tracks
    }    // For all azimuthal angles

    for (std::size_t i = 0; i < nfsrs_; i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double Vi = fsrs_[i]->volume();
      const double Et = mat.Etr(g);
      sflux(g, i) *= 1. / (Vi * Et);
      sflux(g, i) += 4. * PI * src(g, i) / Et;
    }
  }  // For all groups
}

double MOCDriver::calc_keff(const xt::xtensor<double, 2>& flux,
                            const xt::xtensor<double, 2>& old_flux) const {
  double num = 0.;
  double denom = 0.;

  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    const double Vr = fsrs_[i]->volume();
    const auto& mat = *fsrs_[i]->xs();
    for (std::uint32_t g = 0; g < ngroups_; g++) {
      const double VvEf = Vr * mat.vEf(g);
      const double flx = flux(g, i);
      const double oflx = old_flux(g, i);

      num += VvEf * flx;
      denom += VvEf * oflx;
    }
  }

  return keff_ * num / denom;
}

void MOCDriver::fill_source(xt::xtensor<double, 2>& scat_src,
                            const xt::xtensor<double, 2>& flux) const {
  const double inv_k = 1. / keff_;
  const double isotropic = 1. / (4. * PI);
  for (std::uint32_t g = 0; g < ngroups_; g++) {
    for (std::size_t i = 0; i < fsrs_.size(); i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double chi_g = mat.chi(g);
      double Qout = 0.;
      double fiss_rate_i = 0;

      for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
        // Sccatter source
        const double flux_gg_i = flux(gg, i);
        const double Es_gg_to_g = mat.Es_tr(gg, g);
        Qout += Es_gg_to_g * flux_gg_i;

        // Fission source
        const double vEf_gg = mat.vEf(gg);
        Qout += inv_k * chi_g * vEf_gg * flux_gg_i;
      }

      scat_src(g, i) = isotropic * Qout;
    }
  }
}

void MOCDriver::generate_azimuthal_quadrature(std::uint32_t n_angles,
                                              double d) {
  spdlog::info("Creating quadrature");
  spdlog::info("Number of azimuthal angles: {}", n_angles);
  spdlog::info("Maximum track spacing: {} cm", d);

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
    spdlog::debug("Angle {:.5E} pi: weight {:.5E}, width {:.5E}, nx {}, ny {}",
                  ai.phi / PI, ai.wgt, ai.d, ai.nx, ai.ny);
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
        std::pair<UniqueFSR, Vector> fsr_r =
            geometry_->get_fsr_r_local(r_end, u);
        auto ti = geometry_->get_tile_index(r_end, u);
        while (fsr_r.first.fsr && ti) {
          const double d = fsr_r.first.fsr->distance(fsr_r.second, u);
          segments.emplace_back(fsr_r.first.fsr, d,
                                this->get_fsr_indx(fsr_r.first));
          r_end = r_end + d * u;
          ti = geometry_->get_tile_index(r_end, u);
          if (ti) fsr_r = geometry_->get_fsr_r_local(r_end, u);
        }

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
        std::pair<UniqueFSR, Vector> fsr_r =
            geometry_->get_fsr_r_local(r_end, u);
        auto ti = geometry_->get_tile_index(r_end, u);
        while (fsr_r.first.fsr && ti) {
          const double d = fsr_r.first.fsr->distance(fsr_r.second, u);
          segments.emplace_back(fsr_r.first.fsr, d,
                                this->get_fsr_indx(fsr_r.first));
          r_end = r_end + d * u;
          ti = geometry_->get_tile_index(r_end, u);
          if (ti) fsr_r = geometry_->get_fsr_r_local(r_end, u);
        }

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

      // if (tracks.at(i).exit_pos() != comp_tracks.at(ai.ny + i).exit_pos()) {
      //   std::stringstream mssg;
      //   mssg << "Disagreement in track end alignments: "
      //        << tracks.at(i).exit_pos() << " and "
      //        << comp_tracks.at(ai.ny + i).exit_pos() << ".";
      //   spdlog::error(mssg.str());
      //   throw ScarabeeException(mssg.str());
      // }

      tracks.at(i).exit_bc() = this->y_max_bc_;
      comp_tracks.at(ai.nx - 1 - i).exit_bc() = this->y_max_bc_;
    }

    // Go through intersections on bottom side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(ai.ny + i).set_entry_track_flux(
          &comp_tracks.at(i).entry_flux());
      comp_tracks.at(i).set_entry_track_flux(
          &tracks.at(ai.ny + i).entry_flux());

      // if (tracks.at(ai.ny + i).entry_pos() != comp_tracks.at(i).entry_pos())
      // {
      //   std::stringstream mssg;
      //   mssg << "Disagreement in track end alignments: "
      //        << tracks.at(ai.ny + i).entry_pos() << " and "
      //        << comp_tracks.at(i).entry_pos() << ".";
      //   spdlog::error(mssg.str());
      //   throw ScarabeeException(mssg.str());
      // }

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

      // if (tracks.at(i).entry_pos() !=
      //     comp_tracks.at(ai.ny - 1 - i).exit_pos()) {
      //   std::stringstream mssg;
      //   mssg << "Disagreement in track end alignments: "
      //        << tracks.at(i).entry_pos() << " and "
      //        << comp_tracks.at(ai.ny - 1 - i).exit_pos() << ".";
      //   spdlog::error(mssg.str());
      //   throw ScarabeeException(mssg.str());
      // }

      tracks.at(i).entry_bc() = this->x_min_bc_;
      comp_tracks.at(ai.nx + i).exit_bc() = this->x_min_bc_;

      // Right
      tracks.at(ai.nx + i).set_exit_track_flux(
          &comp_tracks.at(nt - 1 - i).entry_flux());
      comp_tracks.at(nt - 1 - i)
          .set_entry_track_flux(&tracks.at(ai.nx + i).exit_flux());

      // if (tracks.at(ai.nx + i).exit_pos() !=
      //     comp_tracks.at(nt - 1 - i).entry_pos()) {
      //   std::stringstream mssg;
      //   mssg << "Disagreement in track end alignments: "
      //        << tracks.at(ai.nx + i).exit_pos() << " and "
      //        << comp_tracks.at(nt - 1 - i).entry_pos() << ".";
      //   spdlog::error(mssg.str());
      //   throw ScarabeeException(mssg.str());
      // }

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
  std::vector<double> approx_vols(nfsrs_, 0.);

  // Go through all angles
  for (std::size_t a = 0; a < angle_info_.size(); a++) {
    const double d = angle_info_[a].d;  // Track width
    auto& tracks = tracks_[a];          // Vector of tracks

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

double MOCDriver::flux(const Vector& r, const Direction& u,
                       std::size_t g) const {
  if (g >= ngroups()) {
    std::stringstream mssg;
    mssg << "Group index g = " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  try {
    const auto& fsr = this->get_fsr(r, u);
    return flux_(g, get_fsr_indx(fsr));
  } catch (ScarabeeException& err) {
    std::stringstream mssg;
    mssg << "Could not find flat source region at r = " << r << " u = " << u
         << "n";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // SHOULD NEVER GET HERE
  auto mssg = "How did I get here ?";
  spdlog::error(mssg);
  throw ScarabeeException(mssg);
  return 0.;
}

double MOCDriver::flux(std::size_t i, std::size_t g) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "FSR index i=" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return flux_(g, i);
}

double MOCDriver::volume(const Vector& r, const Direction& u) const {
  try {
    const auto& fsr = this->get_fsr(r, u);
    return fsr.fsr->volume();
  } catch (ScarabeeException& err) {
    std::stringstream mssg;
    mssg << "Could not find flat source region at r = " << r << " u = " << u
         << "n";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // SHOULD NEVER GET HERE
  auto mssg = "How did I get here ?";
  spdlog::error(mssg);
  throw ScarabeeException(mssg);
  return 0.;
}

double MOCDriver::volume(std::size_t i) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "FSR index i=" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return fsrs_[i]->volume();
}

std::shared_ptr<CrossSection> MOCDriver::xs(const Vector& r,
                                            const Direction& u) const {
  try {
    const auto& fsr = this->get_fsr(r, u);
    return fsr.fsr->xs();
  } catch (ScarabeeException& err) {
    std::stringstream mssg;
    mssg << "Could not find flat source region at r = " << r << " u = " << u
         << "n";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // SHOULD NEVER GET HERE
  auto mssg = "How did I get here ?";
  spdlog::error(mssg);
  throw ScarabeeException(mssg);

  return nullptr;
}

std::shared_ptr<CrossSection> MOCDriver::xs(std::size_t i) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "FSR index i=" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return fsrs_[i]->xs();
}

UniqueFSR MOCDriver::get_fsr(const Vector& r, const Direction& u) const {
  try {
    return geometry_->get_fsr(r, u);
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not find FSR in Cartesian2D.");
    throw err;
  }
}

std::size_t MOCDriver::get_fsr_indx(const UniqueFSR& fsr) const {
  const std::size_t i = fsr_offsets_.at(fsr.fsr->id()) + fsr.instance;
  if (i >= nfsrs_) {
    auto mssg = "FSR index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  return i;
}

void MOCDriver::set_extern_src(const Vector& r, const Direction& u,
                               std::size_t g, double src) {
  std::size_t i;

  try {
    i = this->get_fsr_indx(this->get_fsr(r, u));
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not obtain source region index.");
  }

  if (g >= this->ngroups()) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (src < 0.) {
    auto mssg = "Cannot assign negative external source.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  extern_src_(g, i) = src;
}

double MOCDriver::extern_src(const Vector& r, const Direction& u,
                             std::size_t g) const {
  std::size_t i;

  try {
    i = this->get_fsr_indx(this->get_fsr(r, u));
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not obtain source region index.");
  }

  if (g >= this->ngroups()) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return extern_src_(g, i);
}

void MOCDriver::set_extern_src(std::size_t i, std::size_t g, double src) {
  if (i >= nfsrs_) {
    auto mssg = "Source region index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (g >= this->ngroups()) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (src < 0.) {
    auto mssg = "Cannot assign negative external source.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  extern_src_(g, i) = src;
}

double MOCDriver::extern_src(std::size_t i, std::size_t g) const {
  if (i >= nfsrs_) {
    auto mssg = "Source region index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (g >= this->ngroups()) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return extern_src_(g, i);
}

}  // namespace scarabee
