#include <moc/moc_driver.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>
#include <utils/timer.hpp>
#include <utils/math.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>

namespace scarabee {

MOCDriver::MOCDriver(std::shared_ptr<Cartesian2D> geometry,
                     BoundaryCondition xmin, BoundaryCondition xmax,
                     BoundaryCondition ymin, BoundaryCondition ymax,
                     bool anisotropic)
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
      y_max_bc_(ymax),
      anisotropic_(anisotropic) {
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

  allocate_fsr_data();

  ngroups_ = geometry_->ngroups();

  // get the max-legendre-order for given scattering moments
  if (anisotropic_) {
    max_L_ = 0;
    for (std::size_t i = 0; i < nfsrs_; i++) {
      const auto& mat = *fsrs_[i]->xs();
      const std::size_t l = mat.max_legendre_order();
      if (l > max_L_) max_L_ = l;
    }

    N_lj_ = (max_L_ + 1) * (max_L_ + 1);
  }

  // Allocate arrays and assign indices
  flux_.resize({ngroups_, nfsrs_, N_lj_});
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
  // take all polar angles in case of anisotropic scattering
  if (anisotropic_ == true) {
    n_pol_angles_ *= 2;
  }

  if (n_angles < 4) {
    auto mssg = "MOCDriver must have at least 4 angles.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (n_angles % 4 != 0) {
    // If the number of angles isn't a multiple of 4, we won't be able to make
    // all the boundary condition connections due to an odd number.
    auto mssg = "MOCDriver number of angles must be a multiple of 4.";
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
  trace_tracks();
  segment_renormalization();

  if ((x_min_bc_ == BoundaryCondition::Periodic &&
       x_max_bc_ != BoundaryCondition::Periodic) ||
      (x_min_bc_ != BoundaryCondition::Periodic &&
       x_max_bc_ == BoundaryCondition::Periodic)) {
    auto mssg = "Only one x boundary has a periodic boundary condition.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if ((y_min_bc_ == BoundaryCondition::Periodic &&
       y_max_bc_ != BoundaryCondition::Periodic) ||
      (y_min_bc_ != BoundaryCondition::Periodic &&
       y_max_bc_ == BoundaryCondition::Periodic)) {
    auto mssg = "Only one y boundary has a periodic boundary condition.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  spdlog::info("Determining track connections");
  set_bcs();

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

  if (anisotropic_ == false) {
    // isotropic
    solve_isotropic();
  } else {
    // anisotropic
    solve_anisotropic();
  }

  solved_ = true;

  sim_timer.stop();
  spdlog::info("");
  spdlog::info("Simulation Time: {:.5E} s", sim_timer.elapsed_time());
}

// solve for the isotropic
void MOCDriver::solve_isotropic() {
  flux_.resize({ngroups_, nfsrs_, 1});
  flux_.fill(0.);
  xt::xtensor<double, 2> src;
  src.resize({ngroups_, nfsrs_});
  src.fill(0.);

  // Initialize stabalization matrix (see [1])
  xt::xtensor<double, 2> D;
  D.resize({ngroups_, nfsrs_});
  D.fill(0.);
  for (std::size_t i = 0; i < nfsrs_; i++) {
    const auto xs = this->xs(i);
    for (std::size_t g = 0; g < ngroups_; g++) {
      const double Estr_g_g = xs->Es_tr(g, g);
      if (Estr_g_g < 0.) {
        D(g, i) = -Estr_g_g / xs->Etr(g);
      }
    }
  }

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

    fill_source(src, flux_);
    src += extern_src_;

    // Check for negative source values at beginning of simulation
    bool set_neg_src_to_zero = false;
    if (iteration <= 20) {
      for (std::size_t i = 0; i < src.size(); i++) {
        if (src.flat(i) < 0.) {
          src.flat(i) = 0.;
          set_neg_src_to_zero = true;
        }
      }
    }

    next_flux.fill(0.);
    sweep(next_flux, src);

    // Apply stabalization (see [1])
    for (std::size_t g = 0; g < ngroups_; g++) {
      for (std::size_t i = 0; i < nfsrs_; i++) {
        if (D(g, i) != 0.) {
          next_flux(g, i, 0) += flux_(g, i, 0) * D(g, i);
          next_flux(g, i, 0) /= (1. + D(g, i));
        }
      }
    }

    if (mode_ == SimulationMode::Keff) {
      prev_keff = keff_;
      keff_ = calc_keff(next_flux, flux_);
      rel_diff_keff = std::abs(keff_ - prev_keff) / keff_;
    }

    // Get difference
    max_flx_diff = xt::amax(xt::abs(next_flux - flux_) / next_flux)();

    // Make sure that the flux is positive everywhere !
    bool set_neg_flux_to_zero = false;
    for (std::size_t i = 0; i < next_flux.size(); i++) {
      if (next_flux.flat(i) < 0.) {
        next_flux.flat(i) = 0.;
        set_neg_flux_to_zero = true;
      }
    }

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
    spdlog::info("     iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());

    // Write warnings about negative flux and source
    if (set_neg_src_to_zero) {
      spdlog::warn("Negative source values set to zero");
    }
    if (set_neg_flux_to_zero) {
      spdlog::warn("Negative flux values set to zero");
    }
  }
}

// solve for anisotropic
void MOCDriver::solve_anisotropic() {
  N_lj_ = (max_L_ + 1) * (max_L_ + 1);
  flux_.resize({ngroups_, nfsrs_, N_lj_});
  flux_.fill(0.);

  xt::xtensor<double, 3> src;
  src.resize({ngroups_, nfsrs_, N_lj_});
  src.fill(0.);

  // get the polar angles and azimuthal angle
  // to pre-caluculate the spherical harmonics
  std::span<const double> pol_angl = polar_quad_.polar_angle();
  std::vector<double> polar_angles(pol_angl.begin(), pol_angl.end());
  std::vector<double> azimuthal_angles;
  const std::size_t n_track_angles_ = angle_info_.size();
  azimuthal_angles.reserve(n_track_angles_);
  for (auto& ai : angle_info_) {
    azimuthal_angles.push_back(ai.phi);
  }

  sph_harm_ = SphericalHarmonics(max_L_, azimuthal_angles, polar_angles);

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
      track.entry_flux().fill(1. / std::sqrt(4. * PI));
      track.exit_flux().fill(1. / std::sqrt(4. * PI));
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

    fill_source_anisotropic(src, flux_);
    xt::view(src, xt::all(), xt::all(), 0) += extern_src_;

    // Check for negative zero moment source values at beginning of simulation
    bool set_neg_src_to_zero = false;
    if (iteration <= 20) {
      for (std::size_t g = 0; g < ngroups_; g++) {
        for (std::size_t i = 0; i < nfsrs_; i++) {
          if (src(g, i, 0) < 0.) {
            src(g, i, 0) = 0;
            set_neg_src_to_zero = true;
          }
        }
      }
    }

    next_flux.fill(0.);
    sweep_anisotropic(next_flux, src);

    if (mode_ == SimulationMode::Keff) {
      prev_keff = keff_;
      keff_ = calc_keff(next_flux, flux_);
      rel_diff_keff = std::abs(keff_ - prev_keff) / keff_;
    }

    // Get difference
    auto new_flux = xt::view(next_flux, xt::all(), xt::all(), 0);
    auto old_flux = xt::view(flux_, xt::all(), xt::all(), 0);
    max_flx_diff = xt::amax(xt::abs(new_flux - old_flux) / new_flux)();

    // Make sure that the zero-moment flux is positive everywhere !
    bool set_neg_flux_to_zero = false;
    for (std::size_t g = 0; g < ngroups_; g++) {
      for (std::size_t i = 0; i < nfsrs_; i++) {
        if (next_flux(g, i, 0) < 0.) {
          next_flux(g, i, 0) = 0;
          set_neg_flux_to_zero = true;
        }
      }
    }

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
    spdlog::info("     iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());

    // Write warnings about negative flux and source
    if (set_neg_src_to_zero) {
      spdlog::warn("Negative zero-mometn-source values set to zero.");
    }
    if (set_neg_flux_to_zero) {
      spdlog::warn("Negative zero-moment-flux values set to zero.");
    }
  }
}

void MOCDriver::sweep(xt::xtensor<double, 3>& sflux,
                      const xt::xtensor<double, 2>& src) {
#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(ngroups_); ig++) {
    std::size_t g = static_cast<std::size_t>(ig);

    for (auto& tracks : tracks_) {
      for (std::size_t t = 0; t < tracks.size(); t++) {
        auto& track = tracks[t];
        const double tw = 4. * PI * track.wgt() *
                          track.width();  // Azimuthal weight * track width

        // Get the azimuthal angle (phi) and its cosine for CMFD current
        const Direction u_forw = track.dir();
        const Direction u_back = -u_forw;

        // Load the angular flux for forward direction
        htl::static_vector<double, 6> angflux;
        for (std::size_t p = 0; p < n_pol_angles_; p++)
          angflux.push_back(track.entry_flux()(g, p));

        // Follow track in forward direction
        for (auto& seg : track) {
          const std::size_t i = seg.fsr_indx();
          const double l = seg.length();
          const double Et = seg.xs()->Etr(g);
          const double lEt = l * Et;
          const double Q = src(g, i);
          double delta_sum = 0.;
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            double exp_m1 = mexp(lEt * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux[p] - (Q / Et)) * exp_m1;
            angflux[p] -= delta_flx;
            delta_sum += polar_quad_.wsin()[p] * delta_flx;
          }  // For all polar angles
          sflux(g, i, 0) += tw * delta_sum;
        }  // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.exit_bc() == BoundaryCondition::Vacuum) {
          xt::view(track.exit_track_flux(), g, xt::all()).fill(0.);
        } else {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            track.exit_track_flux()(g, p) = angflux[p];
          }
        }

        // Follow track in backwards direction
        // First, load the backwards angular flux
        for (std::size_t p = 0; p < n_pol_angles_; p++)
          angflux[p] = track.exit_flux()(g, p);

        // Iterate over segments in backwards direction
        for (auto seg_it = track.rbegin(); seg_it != track.rend(); seg_it++) {
          auto& seg = *seg_it;
          const std::size_t i = seg.fsr_indx();
          const double l = seg.length();
          const double Et = seg.xs()->Etr(g);
          const double lEt = l * Et;
          const double Q = src(g, i);
          double delta_sum = 0.;
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            double exp_m1 = mexp(lEt * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux[p] - (Q / Et)) * exp_m1;
            angflux[p] -= delta_flx;
            delta_sum += polar_quad_.wsin()[p] * delta_flx;
          }  // For all polar angles
          sflux(g, i, 0) += tw * delta_sum;
        }  // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.entry_bc() == BoundaryCondition::Vacuum) {
          xt::view(track.entry_track_flux(), g, xt::all()).fill(0.);
        } else {
          for (std::size_t p = 0; p < n_pol_angles_; p++) {
            track.entry_track_flux()(g, p) = angflux[p];
          }
        }
      }  // For all tracks
    }  // For all azimuthal angles

    for (std::size_t i = 0; i < nfsrs_; i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double Vi = fsrs_[i]->volume();
      const double Et = mat.Etr(g);
      sflux(g, i, 0) *= 1. / (Vi * Et);
      sflux(g, i, 0) += 4. * PI * src(g, i) / Et;
    }
  }  // For all groups
}

// anisotropic sweep
void MOCDriver::sweep_anisotropic(xt::xtensor<double, 3>& sflux,
                                  const xt::xtensor<double, 3>& src) {
#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(ngroups_); ig++) {
    std::size_t g = static_cast<std::size_t>(ig);
    for (auto& tracks : tracks_) {
      for (std::size_t t = 0; t < tracks.size(); t++) {
        auto& track = tracks[t];
        htl::static_vector<double, 12> angflux;
        for (std::size_t pp = 0; pp < n_pol_angles_; pp++)
          angflux.push_back(track.entry_flux()(g, pp));
        const double tw = 4. * PI * track.wgt() *
                          track.width();  // Azimuthal weight * track width

        // Follow track in forward direction
        for (auto& seg : track) {
          const std::size_t i = seg.fsr_indx();
          const double l = seg.length();
          const double Et = seg.xs()->Et(g);
          const double lEt = l * Et;
          const std::size_t phi_forward_index = track.phi_index_forward();

          // loop over all polar angles
          std::size_t p = 0;  // index for polar angle
          for (std::size_t pp = 0; pp < n_pol_angles_; pp++) {
            if (pp < n_pol_angles_ / 2) {
              p = pp;
            } else {
              p = pp - n_pol_angles_ / 2;
            }

            double Q = 0.;
            std::span<const double> Y_ljs =
                sph_harm_.spherical_harmonics(phi_forward_index, pp);
            for (std::size_t it_lj = 0; it_lj < N_lj_; it_lj++) {
              Q += src(g, i, it_lj) * Y_ljs[it_lj];
            }

            const double exp_m1 = mexp(lEt * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux[pp] - (Q / Et)) * exp_m1;
            angflux[pp] -= delta_flx;
            const double delta_sum = polar_quad_.wsin()[p] * delta_flx;

            for (std::size_t it_lj = 0; it_lj < N_lj_; it_lj++) {
              sflux(g, i, it_lj) += tw *
                                    (delta_sum + l * Q * polar_quad_.wgt()[p]) *
                                    Y_ljs[it_lj] * 0.5;
            }
          }  // For all polar angles
        }  // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.exit_bc() == BoundaryCondition::Vacuum) {
          xt::view(track.exit_track_flux(), g, xt::all()).fill(0.);
        } else {
          for (std::size_t pp = 0; pp < n_pol_angles_; pp++) {
            track.exit_track_flux()(g, pp) = angflux[pp];
          }
        }

        // Follow track in backwards direction
        for (std::size_t pp = 0; pp < n_pol_angles_; pp++)
          angflux[pp] = track.exit_flux()(g, pp);

        for (auto seg_it = track.rbegin(); seg_it != track.rend(); seg_it++) {
          auto& seg = *seg_it;
          const std::size_t i = seg.fsr_indx();
          const double l = seg.length();
          const double Et = seg.xs()->Et(g);
          const double lEt = l * Et;
          const std::size_t phi_backward_index = track.phi_index_backward();

          // loop over all polar angles
          std::size_t p = 0;
          for (std::size_t pp = 0; pp < n_pol_angles_; pp++) {
            if (pp < n_pol_angles_ / 2) {
              p = pp;
            } else {
              p = pp - n_pol_angles_ / 2;
            }

            // source term evaluation for given azimuthal and polar angle
            double Q = 0.;
            std::span<const double> Y_ljs =
                sph_harm_.spherical_harmonics(phi_backward_index, pp);
            for (std::size_t it_lj = 0; it_lj < N_lj_; it_lj++) {
              Q += src(g, i, it_lj) * Y_ljs[it_lj];
            }

            const double exp_m1 = mexp(lEt * polar_quad_.invs_sin()[p]);
            const double delta_flx = (angflux[pp] - (Q / Et)) * exp_m1;
            angflux[pp] -= delta_flx;
            const double delta_sum = polar_quad_.wsin()[p] * delta_flx;

            for (std::size_t it_lj = 0; it_lj < N_lj_; it_lj++) {
              sflux(g, i, it_lj) += tw *
                                    (delta_sum + l * Q * polar_quad_.wgt()[p]) *
                                    Y_ljs[it_lj] * 0.5;
            }
          }  // For all polar angles
        }  // For all segments along forward direction of track

        // Set incoming flux for next track
        if (track.entry_bc() == BoundaryCondition::Vacuum) {
          xt::view(track.entry_track_flux(), g, xt::all()).fill(0.);
        } else {
          for (std::size_t pp = 0; pp < n_pol_angles_; pp++) {
            track.entry_track_flux()(g, pp) = angflux[pp];
          }
        }
      }  // For all tracks
    }  // For all azimuthal angles

    for (std::size_t i = 0; i < nfsrs_; i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double Vi = fsrs_[i]->volume();
      const double Et = mat.Et(g);
      for (std::size_t it_lj = 0; it_lj < N_lj_; it_lj++) {
        sflux(g, i, it_lj) *= 1. / (Vi * Et);
      }
    }
  }  // For all groups
}

double MOCDriver::calc_keff(const xt::xtensor<double, 3>& flux,
                            const xt::xtensor<double, 3>& old_flux) const {
  double num = 0.;
  double denom = 0.;

#pragma omp parallel
  {
    double num_thrd = 0.;
    double denom_thrd = 0.;

#pragma omp for
    for (int ii = 0; ii < static_cast<int>(fsrs_.size()); ii++) {
      const std::size_t i = static_cast<std::size_t>(ii);
      const double Vr = fsrs_[i]->volume();
      const auto& mat = *fsrs_[i]->xs();
      for (std::uint32_t g = 0; g < ngroups_; g++) {
        const double VvEf = Vr * mat.vEf(g);
        const double flx = flux(g, i, 0);
        const double oflx = old_flux(g, i, 0);

        num_thrd += VvEf * flx;
        denom_thrd += VvEf * oflx;
      }
    }

#pragma omp atomic
    num += num_thrd;

#pragma omp atomic
    denom += denom_thrd;
  }

  return keff_ * num / denom;
}

void MOCDriver::fill_source(xt::xtensor<double, 2>& src,
                            const xt::xtensor<double, 3>& flux) const {
  const double inv_k = 1. / keff_;
  const double isotropic = 1. / (4. * PI);

#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(ngroups_); ig++) {
    const std::size_t g = static_cast<std::size_t>(ig);
    for (std::size_t i = 0; i < fsrs_.size(); i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double chi_g = mat.chi(g);
      double Qout = 0.;

      for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
        // Sccatter source
        const double flux_gg_i = flux(gg, i, 0);
        const double Es_gg_to_g = mat.Es_tr(gg, g);
        Qout += Es_gg_to_g * flux_gg_i;

        // Fission source
        const double vEf_gg = mat.vEf(gg);
        Qout += inv_k * chi_g * vEf_gg * flux_gg_i;
      }

      src(g, i) = isotropic * Qout;
    }
  }
}

void MOCDriver::fill_source_anisotropic(
    xt::xtensor<double, 3>& src, const xt::xtensor<double, 3>& flux) const {
  const double inv_k = 1. / keff_;

#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(ngroups_); ig++) {
    const std::size_t g = static_cast<std::size_t>(ig);
    for (std::size_t i = 0; i < fsrs_.size(); i++) {
      const auto& mat = *fsrs_[i]->xs();
      const double chi_g = mat.chi(g);

      std::size_t it_lj = 0;
      for (std::size_t l = 0; l <= max_L_; l++) {
        for (int j = -static_cast<int>(l); j <= static_cast<int>(l); j++) {
          double Qout = 0.;
          for (std::uint32_t gg = 0; gg < ngroups_; gg++) {
            // Sccatter source
            const double flux_gg_i = flux(gg, i, it_lj);
            const double Es_gg_to_g = mat.Es(l, gg, g);
            Qout += Es_gg_to_g * flux_gg_i;

            // Fission source
            if (l == 0) {
              const double vEf_gg = mat.vEf(gg);
              Qout += inv_k * chi_g * vEf_gg * flux_gg_i;
            }
          }

          src(g, i, it_lj) = Qout;
          it_lj++;

        }  // -l to l
      }  // all scattering moments L
    }  // all flat souce regions
  }  // all groups
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

  angle_info_.resize(n_track_angles_, {0., 0., 0., 0, 0, 0, 0});
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
    angle_info_[i].forward_index = i;
    angle_info_[i].backward_index = i + n_track_angles_;
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

void MOCDriver::trace_tracks() {
  spdlog::info("Tracing tracks");

  std::uint32_t n_track_angles_ =
      static_cast<std::uint32_t>(angle_info_.size());
  const double Dx = geometry_->x_max() - geometry_->x_min();
  const double Dy = geometry_->y_max() - geometry_->y_min();

  tracks_.resize(n_track_angles_);

#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(n_track_angles_); ii++) {
    const std::size_t i = static_cast<std::size_t>(ii);

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

        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt, ai.d,
                                segments, ai.forward_index, ai.backward_index);

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

        tracks_[i].emplace_back(r_start, r_end, u, ai.phi, ai.wgt, ai.d,
                                segments, ai.forward_index, ai.backward_index);

        if (t < ai.nx) {
          x += dx;
        } else {
          y += dy;
        }
      }
    }
  }
}

std::vector<std::pair<std::size_t, double>> MOCDriver::trace_fsr_segments(
    const Vector r_start, const Direction& u) const {
  std::vector<std::pair<std::size_t, double>> out;

  Vector r_end = r_start;
  std::pair<UniqueFSR, Vector> fsr_r = geometry_->get_fsr_r_local(r_end, u);
  auto ti = geometry_->get_tile_index(r_end, u);
  while (fsr_r.first.fsr && ti) {
    const double d = fsr_r.first.fsr->distance(fsr_r.second, u);
    out.emplace_back(this->get_fsr_indx(fsr_r.first), d);
    r_end = r_end + d * u;
    ti = geometry_->get_tile_index(r_end, u);
    if (ti) fsr_r = geometry_->get_fsr_r_local(r_end, u);
  }

  return out;
}

void MOCDriver::set_periodic_bcs_x() {
  // We are setting the boundarys a x_min and x_max, so we need to do all
  // the y tracks that start and end on those sides.
  for (std::size_t a = 0; a < angle_info_.size(); a++) {
    const auto& ai = angle_info_[a];
    auto& tracks = tracks_[a];

    if (ai.phi < PI_2) {
      for (std::size_t j = 0; j < ai.ny; j++) {
        tracks[ai.nx + j].set_exit_track_flux(&tracks[j].entry_flux());
        tracks[j].set_entry_track_flux(&tracks[ai.nx + j].exit_flux());

        tracks[ai.nx + j].exit_bc() = BoundaryCondition::Periodic;
        tracks[j].entry_bc() = BoundaryCondition::Periodic;
      }
    } else {
      for (std::size_t j = 0; j < ai.ny; j++) {
        tracks[j].set_exit_track_flux(&tracks[ai.nx + j].entry_flux());
        tracks[ai.nx + j].set_entry_track_flux(&tracks[j].exit_flux());

        tracks[j].exit_bc() = BoundaryCondition::Periodic;
        tracks[ai.nx + j].entry_bc() = BoundaryCondition::Periodic;
      }
    }
  }
}

void MOCDriver::set_periodic_bcs_y() {
  // We are setting the boundarys a y_min and y_max, so we need to do all
  // the x tracks that start and end on those sides.
  for (std::size_t a = 0; a < angle_info_.size(); a++) {
    const auto& ai = angle_info_[a];
    auto& tracks = tracks_[a];

    if (ai.phi < PI_2) {
      for (std::size_t i = 0; i < ai.nx; i++) {
        tracks[i].set_exit_track_flux(&tracks[ai.ny + i].entry_flux());
        tracks[ai.ny + i].set_entry_track_flux(&tracks[i].exit_flux());

        tracks[i].exit_bc() = BoundaryCondition::Periodic;
        tracks[ai.ny + i].entry_bc() = BoundaryCondition::Periodic;
      }
    } else {
      for (std::size_t i = 0; i < ai.nx; i++) {
        tracks[i].set_entry_track_flux(&tracks[ai.ny + i].exit_flux());
        tracks[ai.ny + i].set_exit_track_flux(&tracks[i].entry_flux());

        tracks[i].entry_bc() = BoundaryCondition::Periodic;
        tracks[ai.ny + i].exit_bc() = BoundaryCondition::Periodic;
      }
    }
  }
}

void MOCDriver::set_ref_vac_bcs_y_max() {
  // We only go through the first half of the tracks, where phi < pi / 2.
  for (std::size_t a = 0; a < angle_info_.size() / 2; a++) {
    const auto& ai = angle_info_[a];
    auto& tracks = tracks_[a];
    auto& comp_tracks = *(tracks_.rbegin() + a);

    // Go through intersections on top side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(i).set_exit_track_flux(&comp_tracks.at(ai.ny + i).exit_flux());
      comp_tracks.at(ai.ny + i).set_exit_track_flux(&tracks.at(i).exit_flux());

      tracks.at(i).exit_bc() = this->y_max_bc_;
      comp_tracks.at(ai.ny + i).exit_bc() = this->y_max_bc_;
    }
  }
}

void MOCDriver::set_ref_vac_bcs_y_min() {
  // We only go through the first half of the tracks, where phi < pi / 2.
  for (std::size_t a = 0; a < angle_info_.size() / 2; a++) {
    const auto& ai = angle_info_[a];
    auto& tracks = tracks_[a];
    auto& comp_tracks = *(tracks_.rbegin() + a);

    // Go through intersections on bottom side
    for (std::uint32_t i = 0; i < ai.nx; i++) {
      tracks.at(ai.ny + i).set_entry_track_flux(
          &comp_tracks.at(i).entry_flux());
      comp_tracks.at(i).set_entry_track_flux(
          &tracks.at(ai.ny + i).entry_flux());

      tracks.at(ai.ny + i).entry_bc() = this->y_min_bc_;
      comp_tracks.at(i).entry_bc() = this->y_min_bc_;
    }
  }
}

void MOCDriver::set_ref_vac_bcs_x_max() {
  // We only go through the first half of the tracks, where phi < pi / 2.
  for (std::size_t a = 0; a < angle_info_.size() / 2; a++) {
    const auto& ai = angle_info_[a];
    const std::uint32_t nt = ai.nx + ai.ny;
    auto& tracks = tracks_[a];
    auto& comp_tracks = *(tracks_.rbegin() + a);

    // Go down right side
    for (std::uint32_t i = 0; i < ai.ny; i++) {
      tracks.at(ai.nx + i).set_exit_track_flux(
          &comp_tracks.at(nt - 1 - i).entry_flux());
      comp_tracks.at(nt - 1 - i)
          .set_entry_track_flux(&tracks.at(ai.nx + i).exit_flux());

      tracks.at(ai.nx + i).exit_bc() = this->x_max_bc_;
      comp_tracks.at(nt - 1 - i).entry_bc() = this->x_max_bc_;
    }
  }
}

void MOCDriver::set_ref_vac_bcs_x_min() {
  // We only go through the first half of the tracks, where phi < pi / 2.
  for (std::size_t a = 0; a < angle_info_.size() / 2; a++) {
    const auto& ai = angle_info_[a];
    auto& tracks = tracks_[a];
    auto& comp_tracks = *(tracks_.rbegin() + a);

    // Go down left side
    // Go down left/right sides
    for (std::uint32_t i = 0; i < ai.ny; i++) {
      tracks.at(i).set_entry_track_flux(
          &comp_tracks.at(ai.ny - 1 - i).exit_flux());
      comp_tracks.at(ai.ny - 1 - i)
          .set_exit_track_flux(&tracks.at(i).entry_flux());

      tracks.at(i).entry_bc() = this->x_min_bc_;
      comp_tracks.at(ai.ny - 1 - i).exit_bc() = this->x_min_bc_;
    }
  }
}

void MOCDriver::set_bcs() {
  if (x_min_bc_ == BoundaryCondition::Periodic) {
    set_periodic_bcs_x();
  } else {
    set_ref_vac_bcs_x_max();
    set_ref_vac_bcs_x_min();
  }

  if (y_min_bc_ == BoundaryCondition::Periodic) {
    set_periodic_bcs_y();
  } else {
    set_ref_vac_bcs_y_max();
    set_ref_vac_bcs_y_min();
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

void MOCDriver::allocate_fsr_data() {
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

double MOCDriver::flux(const Vector& r, const Direction& u, std::size_t g,
                       std::size_t lj) const {
  if (g >= ngroups()) {
    std::stringstream mssg;
    mssg << "Group index g = " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (lj >= N_lj_) {
    std::stringstream mssg;
    mssg << "Spherical harmonic index lj = " << lj << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  try {
    const auto& fsr = this->get_fsr(r, u);
    return flux_(g, get_fsr_indx(fsr), lj);
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

double MOCDriver::flux(std::size_t i, std::size_t g, std::size_t lj) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "FSR index i=" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g >= ngroups()) {
    std::stringstream mssg;
    mssg << "Group index g = " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (lj >= N_lj_) {
    std::stringstream mssg;
    mssg << "Spherical harmonic index lj = " << lj << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return flux_(g, i, lj);
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

const std::shared_ptr<CrossSection>& MOCDriver::xs(const Vector& r,
                                                   const Direction& u) const {
  UniqueFSR fsr;
  try {
    fsr = this->get_fsr(r, u);
  } catch (ScarabeeException& err) {
    std::stringstream mssg;
    mssg << "Could not find flat source region at r = " << r << " u = " << u
         << "n";
    err.add_to_exception(mssg.str());
    throw err;
  }

  return fsr.fsr->xs();
}

const std::shared_ptr<CrossSection>& MOCDriver::xs(std::size_t i) const {
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

std::vector<std::size_t> MOCDriver::get_all_fsr_in_cell(
    const Vector& r, const Direction& u) const {
  auto fsrs = geometry_->get_all_fsr_in_cell(r, u);
  std::vector<std::size_t> inds(fsrs.size(), 0);

  for (std::size_t i = 0; i < fsrs.size(); i++) {
    inds[i] = this->get_fsr_indx(fsrs[i]);
  }

  return inds;
}

void MOCDriver::set_extern_src(const Vector& r, const Direction& u,
                               std::size_t g, double src) {
  std::size_t i = 0;

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
  std::size_t i = 0;

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

std::shared_ptr<CrossSection> MOCDriver::homogenize() const {
  const std::size_t NR = this->nfsr();
  std::vector<std::size_t> regions(NR, 0);
  for (std::size_t i = 0; i < NR; i++) {
    regions[i] = i;
  }

  return this->homogenize(regions);
}

std::shared_ptr<CrossSection> MOCDriver::homogenize(
    const std::vector<std::size_t>& regions) const {
  // Check all regions are valid
  if (regions.size() > this->nregions()) {
    auto mssg =
        "The number of provided regions is greater than the number of regions.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto m : regions) {
    if (m >= this->nfsr()) {
      auto mssg = "Invalid region index in homogenization list.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // We now begin homogenization
  const std::size_t NR = regions.size();
  const std::size_t NG = this->ngroups();

  std::size_t max_l = 0;
  for (const auto m : regions) {
    const auto m_max_l = this->xs(m)->max_legendre_order();
    if (m_max_l > max_l) {
      max_l = m_max_l;
    }
  }

  xt::xtensor<double, 1> Et = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Dtr = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 3> Es = xt::zeros<double>({max_l + 1, NG, NG});
  xt::xtensor<double, 1> Ef = xt::zeros<double>({NG});
  xt::xtensor<double, 1> vEf = xt::zeros<double>({NG});
  xt::xtensor<double, 1> chi = xt::zeros<double>({NG});

  // We need to calculate the total fission production in each volume for
  // generating the homogenized fission spectrum.
  std::vector<double> fiss_prod(NR, 0.);
  std::size_t j = 0;
  for (const auto i : regions) {
    const auto& mat = this->xs(i);
    const double V = this->volume(i);
    for (std::size_t g = 0; g < NG; g++) {
      fiss_prod[j] += mat->vEf(g) * flux(i, g) * V;
    }
    j++;
  }
  const double sum_fiss_prod =
      std::accumulate(fiss_prod.begin(), fiss_prod.end(), 0.);
  const double invs_sum_fiss_prod =
      sum_fiss_prod > 0. ? 1. / sum_fiss_prod : 1.;

  for (std::size_t g = 0; g < NG; g++) {
    // Get the sum of flux*volume for this group
    double sum_fluxV = 0.;
    for (const auto i : regions) {
      sum_fluxV += this->flux(i, g) * this->volume(i);
    }
    const double invs_sum_fluxV = 1. / sum_fluxV;

    j = 0;
    for (const auto i : regions) {
      const auto& mat = this->xs(i);
      const double V = volume(i);
      const double flx = flux(i, g);
      const double coeff = invs_sum_fluxV * flx * V;
      Dtr(g) += coeff * mat->Dtr(g);
      Ea(g) += coeff * mat->Ea(g);
      Ef(g) += coeff * mat->Ef(g);
      vEf(g) += coeff * mat->vEf(g);

      chi(g) += invs_sum_fiss_prod * fiss_prod[j] * mat->chi(g);

      for (std::size_t l = 0; l <= max_l; l++) {
        for (std::size_t gg = 0; gg < NG; gg++) {
          Es(l, g, gg) += coeff * mat->Es(l, g, gg);
        }
      }

      j++;
    }

    // Reconstruct total xs from absorption and scattering
    Et(g) = Ea(g) + xt::sum(xt::view(Es, 0, g, xt::all()))();
  }

  return std::make_shared<CrossSection>(Et, Dtr, Ea, Es, Ef, vEf, chi);
}

xt::xtensor<double, 1> MOCDriver::homogenize_flux_spectrum() const {
  const std::size_t NR = this->nfsr();
  std::vector<std::size_t> regions(NR, 0);
  for (std::size_t i = 0; i < NR; i++) {
    regions[i] = i;
  }

  return this->homogenize_flux_spectrum(regions);
}

xt::xtensor<double, 1> MOCDriver::homogenize_flux_spectrum(
    const std::vector<std::size_t>& regions) const {
  // Check all regions are valid
  if (regions.size() > this->nfsr()) {
    auto mssg =
        "The number of provided regions is greater than the number of regions.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto m : regions) {
    if (m >= this->nregions()) {
      auto mssg = "Invalid region index in homogenization list.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  const std::size_t NG = this->ngroups();

  // First, calculate the sum of the volumes
  double sum_V = 0.;
  for (const auto i : regions) {
    sum_V += this->volume(i);
  }
  const double invs_sum_V = 1. / sum_V;

  xt::xtensor<double, 1> spectrum = xt::zeros<double>({NG});
  for (std::size_t g = 0; g < NG; g++) {
    for (const auto i : regions) {
      spectrum(g) += invs_sum_V * this->volume(i) * this->flux(i, g);
    }
  }

  return spectrum;
}

void MOCDriver::apply_criticality_spectrum(const xt::xtensor<double, 1>& flux) {
  if (this->solved() == false) {
    auto mssg =
        "Cannot apply criticality spectrum when problem has not been solved.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (flux.size() != this->ngroups()) {
    auto mssg =
        "Length of criticality flux spectrum does not agree with the number of "
        "groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  xt::xtensor<double, 1> group_mult = this->homogenize_flux_spectrum();
  group_mult = 1. / group_mult;
  group_mult *= flux;

  for (std::size_t g = 0; g < this->ngroups(); g++) {
    for (std::size_t i = 0; i < this->nfsr(); i++) {
      flux_(g, i, 0) *= group_mult(g);
    }
  }
}

void MOCDriver::save_bin(const std::string& fname) const {
  if (std::filesystem::exists(fname)) {
    std::filesystem::remove(fname);
  }

  std::ofstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryOutputArchive arc(file);

  arc(*this);
}

std::shared_ptr<MOCDriver> MOCDriver::load_bin(const std::string& fname) {
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  std::shared_ptr<MOCDriver> out(new MOCDriver());

  std::ifstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryInputArchive arc(file);

  arc(*out);

  return out;
}

}  // namespace scarabee

// REFERENCES
// [1] G. Gunow, B. Forget, and K. Smith, “Stabilization of multi-group neutron
//     transport with transport-corrected cross-sections,” Ann. Nucl. Energy,
//     vol. 126, pp. 211–219, 2019, doi: 10.1016/j.anucene.2018.10.036.