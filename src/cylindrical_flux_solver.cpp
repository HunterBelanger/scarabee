#include <cylindrical_flux_solver.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <xtensor/xbuilder.hpp>

namespace scarabee {

CylindricalFluxSolver::CylindricalFluxSolver(
    std::shared_ptr<CylindricalCell> cell)
    : flux_(),
      j_ext_(),
      x_(),
      cell_(cell),
      k_(1.),
      a_(1.),
      k_tol_(1.E-5),
      flux_tol_(1.E-5),
      solved_(false) {
  if (cell_ == nullptr) {
    auto mssg = "Provided CylindricalCell was a nullptr.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Initialize flux with 1's everywhere
  flux_.resize({ngroups(), nregions()});
  flux_.fill(1.);
  j_ext_.resize({ngroups()});
  j_ext_.fill(0.);
  x_.resize({ngroups()});
  x_.fill(0.);
}

void CylindricalFluxSolver::set_albedo(double a) {
  if (a < 0. || a > 1.) {
    auto mssg = "Albedo must be in the interval [0,1].";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  a_ = a;
  solved_ = false;
}

void CylindricalFluxSolver::set_flux_tolerance(double ftol) {
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
  solved_ = false;
}

void CylindricalFluxSolver::set_keff_tolerance(double ktol) {
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

  k_tol_ = ktol;
  solved_ = false;
}

double CylindricalFluxSolver::Qscat(std::uint32_t g, std::size_t i,
                                    const xt::xtensor<double, 2>& flux) const {
  double Qout = 0.;
  const auto& xs = cell_->xs(i);

  for (std::uint32_t gg = 0; gg < ngroups(); gg++) {
    const double flux_gg_i = flux(gg, i);

    // Scattering into group g, excluding g -> g
    if (gg != g) Qout += xs->Es_tr(gg, g) * flux_gg_i;
  }

  return Qout;
}

double CylindricalFluxSolver::Qfiss(std::uint32_t g, std::size_t i,
                                    const xt::xtensor<double, 2>& flux) const {
  double Qout = 0.;
  const double inv_k = 1. / k_;
  const auto& xs = cell_->xs(i);
  const double chi_g = xs->chi(g);

  for (std::uint32_t gg = 0; gg < ngroups(); gg++) {
    Qout += xs->vEf(gg) * flux(gg, i);
  }

  Qout *= inv_k * chi_g;

  return Qout;
}

double CylindricalFluxSolver::calc_keff(
    const xt::xtensor<double, 2>& flux) const {
  double keff = 0.;
  for (std::size_t r = 0; r < nregions(); r++) {
    const double Vr = cell_->volume(r);
    const auto& xs = cell_->xs(r);
    for (std::uint32_t g = 0; g < ngroups(); g++) {
      keff += Vr * xs->vEf(g) * flux(g, r);
    }
  }

  return keff;
}

double CylindricalFluxSolver::calc_flux_rel_diff(
    const xt::xtensor<double, 2>& flux,
    const xt::xtensor<double, 2>& next_flux) const {
  double max_rel_diff = 0.;

  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::size_t r = 0; r < nregions(); r++) {
      const double rel_diff =
          std::abs(next_flux(g, r) - flux(g, r)) / flux(g, r);

      if (rel_diff > max_rel_diff) max_rel_diff = rel_diff;
    }
  }

  return max_rel_diff;
}

void CylindricalFluxSolver::fill_fission_source(
    xt::xtensor<double, 2>& source, const xt::xtensor<double, 2>& flux) const {
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::size_t r = 0; r < nregions(); r++) {
      source(g, r) = Qfiss(g, r, flux);
    }
  }
}

void CylindricalFluxSolver::fill_scatter_source(
    xt::xtensor<double, 2>& source, const xt::xtensor<double, 2>& flux) const {
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::size_t r = 0; r < nregions(); r++) {
      source(g, r) = Qscat(g, r, flux);
    }
  }
}

void CylindricalFluxSolver::solve() {
  spdlog::info("Solving for keff.");
  spdlog::info("keff tolerance: {:.5E}", k_tol_);
  spdlog::info("Flux tolerance: {:.5E}", flux_tol_);

  // Make sure the cell is solved
  if (cell_->solved() == false) {
    auto mssg = "Cannot solve for flux if cell is not solved.";
    throw ScarabeeException(mssg);
  }

  // Create a new array to hold the source according to the current flux, and
  // another for the next generation flux.
  xt::xtensor<double, 2> scat_source(flux_.shape());
  xt::xtensor<double, 2> fiss_source(flux_.shape());
  xt::xtensor<double, 2> next_flux(flux_.shape());
  k_ = calc_keff(flux_);
  double old_keff = 100.;
  spdlog::info("Initial keff {:.5f}", k_);

  std::size_t outer_iter = 0;
  // Outer Generations
  while (std::abs(old_keff - k_) > k_tol_) {
    outer_iter++;

    // At the begining of a generation, we calculate the fission source
    fill_fission_source(fiss_source, flux_);

    // Copy flux into next_flux, so that we can continually use next_flux
    // in the inner iterations for the scattering source.
    next_flux = flux_;

    double max_flux_diff = 100.;
    std::size_t inner_iter = 0;
    // Inner Iterations
    while (max_flux_diff > flux_tol_) {
      inner_iter++;

      // At the begining of an inner iteration, we calculate the fission source
      fill_scatter_source(scat_source, next_flux);

      // From the sources, we calculate the new flux values
      for (std::uint32_t g = 0; g < ngroups(); g++) {
        for (std::size_t r = 0; r < nregions(); r++) {
          const double Yr = cell_->Y(a_, g, r);

          double Xr = 0.;
          for (std::size_t k = 0; k < nregions(); k++) {
            Xr +=
                (fiss_source(g, k) + scat_source(g, k)) * cell_->X(a_, g, r, k);
          }

          next_flux(g, r) = Xr + j_ext_[g] * Yr;
        }
      }

      // Calculate the max difference in the flux
      max_flux_diff = calc_flux_rel_diff(flux_, next_flux);
      spdlog::debug("Inner iteration {} flux diff {:.5E}", inner_iter,
                    max_flux_diff);

      // Copy next_flux into flux for calculating next relative difference
      flux_ = next_flux;
    }  // End of Inner Iterations

    // Calculate keff
    old_keff = k_;
    k_ = calc_keff(next_flux);
    spdlog::info("Iteration {} keff {:.5f}", outer_iter, k_);

    // Assign next_flux to be the flux
    std::swap(next_flux, flux_);
  }  // End of Outer Generations

  // Now that we have the solution, we need to get the number of source
  // neutrons reaching the boundary, x, from Stamm'ler and Abbate.
  // These are used when calculating the currents.
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::size_t r = 0; r < nregions(); r++) {
      x_[g] += (Qfiss(g, r, flux_) + Qscat(g, r, flux_)) * cell_->x(g, r);
    }
  }

  solved_ = true;
}

}  // namespace scarabee
