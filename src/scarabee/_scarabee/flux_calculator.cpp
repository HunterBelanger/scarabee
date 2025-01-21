#include <data/flux_calculator.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <cmath>

namespace scarabee {

FluxCalculator::FluxCalculator(const xt::xtensor<double, 1>& energy_boundaries,
                               const xt::xtensor<double, 1>& sig_t_r,
                               const xt::xtensor<double, 1>& sig_s_r,
                               double awr_r)
    : energy_boundaries_(energy_boundaries),
      avg_energy_(),
      dlt_energy_(),
      flux_(),
      sig_t_r_(sig_t_r),
      sig_s_r_(sig_s_r),
      chi_(),
      background_nuclides_(),
      A_r_(awr_r),
      a_r_(this->a_from_awr(A_r_)),
      g_max_r(0) {
  if (sig_t_r_.size() != energy_boundaries_.size() - 1) {
    auto mssg =
        "Number of total cross sections does not agree with the number of "
        "energies.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (sig_s_r_.size() != energy_boundaries_.size() - 1) {
    auto mssg =
        "Number of scattering cross sections does not agree with the number of "
        "energies.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Create the additional energy arrays
  avg_energy_ = xt::zeros<double>({sig_t_r_.size()});
  dlt_energy_ = xt::zeros<double>({sig_t_r_.size()});
  for (std::size_t g = 0; g < avg_energy_.size(); g++) {
    double Elow = energy_boundaries_(g);
    double Ehi = energy_boundaries_(g + 1);
    avg_energy_(g) = std::sqrt(Elow * Ehi);
    dlt_energy_(g) = Ehi - Elow;
  }

  // Create the fission spectrum
  generate_watt_spectrum();

  // Initialize empty flux
  flux_ = xt::zeros<double>({sig_t_r_.size()});
}

void FluxCalculator::generate_watt_spectrum() {
  constexpr double a = 988000.0;
  constexpr double b = 2.249E-6;
  constexpr double invs_a = 1. / a;

  chi_ = xt::zeros<double>({avg_energy_.size()});

  double integral = 0.;
  for (std::size_t g = 0; g < chi_.size(); g++) {
    const double E = avg_energy_(g);
    chi_(g) = std::exp(-E * invs_a) * std::sinh(std::sqrt(b * E));
    integral += chi_(g) * dlt_energy_(g);
  }

  chi_ /= integral;
}

void FluxCalculator::add_background_nuclide(double sig_b, double awr) {
  if (sig_b <= 0.) {
    auto mssg = "Nuclide's background cross section must be greater than zero.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (awr <= 0.) {
    auto mssg = "Nuclide's atomic weight ratio must be greater than zero.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  background_nuclides_.push_back({awr, a_from_awr(awr), sig_b, 0});
}

void FluxCalculator::solve() {
  const std::size_t NG = avg_energy_.size();

  // Compute total background xs which is constant in energy
  double sig_b = 0.;
  for (const auto& n : background_nuclides_) sig_b += n.sig_b;

  // Sett the starting g_max indices to be NG-1 for all nuclides
  g_max_r = NG - 1;
  for (auto& n : background_nuclides_) n.g_max = NG - 1;

  // We work from HIGH energy (end of arrays) to LOW energy (beginning of
  // arrays). To start, we set the flux at the highest energy
  flux_(NG - 1) = chi_(NG - 1) / (sig_t_r_(NG - 1) + sig_b);

  // Source in previous group
  double Sgm1 = 0.;

  // Iterate over energy BACKWARDS
  for (int ig = static_cast<int>(NG - 2); ig >= 0; ig--) {
    const std::size_t g = static_cast<std::size_t>(ig);

    // Get total xs at our energy
    const double Et_g = sig_t_r_(g) + sig_b;

    // Compute source
    double Sg = compute_scatter_source(Sgm1, g);

    // Compute flux
    flux_(g) = (Sg + chi_(g)) / Et_g;

    // Set previous source for next iteration
    Sgm1 = Sg;
  }
}

double FluxCalculator::compute_scatter_source(double Sgm1, std::size_t g) {
  const double Eg = avg_energy_(g);

  double Sg = Sgm1;

  // Add contribution for the previous group of all nuclides.
  // Previous group means go up an energy ! Hence the g+1 !
  const double flux_gm1 = flux_(g + 1);
  const double Dlt_Egm1 = dlt_energy_(g + 1);
  const double Egm1 = avg_energy_(g + 1);

  Sg += sig_s_r_(g + 1) * flux_gm1 * Dlt_Egm1 / ((1. - a_r_) * Egm1);
  for (const auto& n : background_nuclides_) {
    Sg += n.sig_b * flux_gm1 * Dlt_Egm1 / ((1. - n.a) * Egm1);
  }

  // We now remove the no-longer needed contributions from energies that
  // are too high to reach our current energy. Start with that from the
  // resonant isotope first.
  double Emax_prev = Egm1 / a_r_;
  double Emax = Eg / a_r_;
  while (avg_energy_(g_max_r) > Emax_prev) g_max_r -= 1;
  for (std::size_t gg = g_max_r; gg > g; gg--) {
    if (avg_energy_(gg) > Emax) {
      Sg -= sig_s_r_(gg) * flux_(gg) * dlt_energy_(gg) /
            ((1. - a_r_) * avg_energy_(gg));
    } else {
      break;
    }
  }

  // Do components from other nuclides
  for (auto& n : background_nuclides_) {
    Emax_prev = Egm1 / n.a;
    Emax = Eg / n.a;
    while (avg_energy_(n.g_max) > Emax_prev) n.g_max -= 1;
    for (std::size_t gg = n.g_max; gg > g; gg--) {
      if (avg_energy_(gg) > Emax) {
        Sg -= n.sig_b * flux_(gg) * dlt_energy_(gg) /
              ((1. - n.a) * avg_energy_(gg));
      } else {
        break;
      }
    }
  }

  return Sg;
}

}  // namespace scarabee
