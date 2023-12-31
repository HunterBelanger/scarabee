#include <cylindrical_flux_solver.hpp>
#include <utils/scarabee_exception.hpp>

CylindricalFluxSolver::CylindricalFluxSolver(
    std::shared_ptr<CylindricalCell> cell)
    : flux_(), j_ext_(), cell_(cell), k_(1.), a_(1.) {
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

double CylindricalFluxSolver::calc_flux_component(std::uint32_t g,
                                                  std::size_t i) const {
  const double Yi = cell_->Y(a_, g, i);
  double Xi = 0.;

  for (std::size_t k = 0; k < nregions(); k++) {
    Xi += Q(g, k) * cell_->X(a_, g, i, k);
  }

  return Xi + j_ext_[g] * Yi;
}
