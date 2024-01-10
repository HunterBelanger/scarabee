#include <cylindrical_flux_solver.hpp>
#include <utils/scarabee_exception.hpp>

#include <iostream>

CylindricalFluxSolver::CylindricalFluxSolver(
    std::shared_ptr<CylindricalCell> cell)
    : flux_(), j_ext_(), cell_(cell), k_(1.), a_(1.), k_tol_(1.E-5), solved_(false) {
  if (cell_ == nullptr) {
    throw ScarabeeException("Provided CylindricalCell was a nullptr.");
  }

  flux_.reallocate({ngroups(), nregions()});
  flux_.fill(1.);  // Initialize flux with 1's everywhere
  j_ext_.resize(ngroups(), 0.);
}

void CylindricalFluxSolver::set_albedo(double a) {
  if (a < 0. || a > 1.) {
    throw ScarabeeException("Albedo must be in the interval [0,1].");
  }

  a_ = a;
}

void CylindricalFluxSolver::set_keff_tolerance(double ktol) {
  if (ktol <= 0.) {
    throw ScarabeeException("Tolerance for keff must be in the interval (0., 0.1).");
  }

  if (ktol >= 0.1) {
    throw ScarabeeException("Tolerance for keff must be in the interval (0., 0.1).");
  }

  k_tol_ = ktol;
}

double CylindricalFluxSolver::Q(std::uint32_t g, std::size_t i) const {
  double Qout = 0.;
  const double inv_k = 1. / k_;
  const auto& mat = cell_->mat(i);
  const double chi_prompt_g = mat.chi_prompt(g);

  for (std::uint32_t gg = 0; gg < ngroups(); gg++) {
    const double Ef_gg = mat.Ef(gg);
    const double flux_gg_i = flux(gg, i);

    // Prompt Fission
    Qout += inv_k * chi_prompt_g * mat.nu_prompt(gg) * Ef_gg * flux_gg_i;

    // Delayed Fission
    for (std::uint32_t f = 0; f < mat.ndelayed_families(); f++) {
      const double chi_delayed_g_f = mat.chi_delayed(g, f);
      Qout +=
          inv_k * chi_delayed_g_f * mat.nu_delayed(gg, f) * Ef_gg * flux_gg_i;
    }

    // Scattering into group g, excluding g -> g
    if (gg != g) Qout += mat.Es_tr(gg, g) * flux_gg_i;
  }

  return Qout;
}

double CylindricalFluxSolver::calc_keff(const NDArray<double>& flux) const {
  double keff = 0.;
  for (std::size_t r = 0; r < nregions(); r++) {
    const double Vr = cell_->V(r);
    const auto& mat = cell_->mat(r);
    for (std::uint32_t g = 0; g < ngroups(); g++) {
      keff += Vr * mat.nu(g) * mat.Ef(g) * flux(g, r);
    }
  }

  return keff;
}

void CylindricalFluxSolver::solve() {
  // Create a new array to hold the source according to the current flux, and
  // another for the next generation flux.
  NDArray<double> source(flux_.shape());
  NDArray<double> next_flux(flux_.shape());
  k_ = calc_keff(flux_);
  //k_ = 1.;
  double old_keff = 100.;
  std::size_t Ngenerations = 0;

  while (std::abs(old_keff - k_) > k_tol_) {
    Ngenerations++;

    // First, we fill the source
    for (std::uint32_t g = 0; g < ngroups(); g++) {
      for (std::size_t r = 0; r < nregions(); r++) {
        source(g, r) = Q(g, r);
      }
    }

    // From the source, we calculate the new flux values
    for (std::uint32_t g = 0; g < ngroups(); g++) {
      for (std::size_t r = 0; r < nregions(); r++) {
        const double Yr = cell_->Y(a_, g, r);
        double Xr = 0.;

        for (std::size_t k = 0; k < nregions(); k++) {
          Xr += source(g, k) * cell_->X(a_, g, r, k);
        }
        next_flux(g, r) = Xr + j_ext_[g] * Yr;
      }
    }

    // Calculate keff
    //double k_num = 0.;
    //double k_denom = 1.;
    //for (std::size_t i = 0; i < flux_.size(); i++) {
    //  k_num += next_flux[i] * flux_[i];
    //  k_denom += flux_[i] * flux_[i];
    //}
    //old_keff = k_;
    //k_ = std::sqrt(k_num / k_denom);
    old_keff = k_;
    k_ = calc_keff(next_flux);

    // Normalize next flux by keff
    for (std::size_t i = 0; i < next_flux.size(); i++) {
      next_flux[i] /= k_;
    }

    // Assign next_flux to be the flux
    std::swap(next_flux.data_vector(), flux_.data_vector());

    std::cout << " Gen: " << Ngenerations << "  keff = " << std::fixed << k_ << "\n";
  }
}
