#include <reflector_sn.hpp>

#include <utils/math.hpp>
#include <utils/logging.hpp>

#include <xtensor/xio.hpp>

namespace scarabee {

void ReflectorSN::solve() {
  const std::size_t NG = xs_[0]->ngroups();
  const std::size_t NR = xs_.size();

  flux_ = xt::ones<double>({NG, NR});
  xt::xtensor<double, 2> inner_flux = flux_;
  xt::xtensor<double, 2> next_inner_flux = xt::zeros<double>({NG, NR});
  xt::xtensor<double, 2> old_outer_flux = xt::zeros<double>({NG, NR});
  xt::xtensor<double, 2> Qfiss = xt::zeros<double>({NG, NR});
  xt::xtensor<double, 2> Qscat = xt::zeros<double>({NG, NR});
  xt::xtensor<double, 2> Q = xt::zeros<double>({NG, NR});

  keff_ = 1.;

  // Outer Iterations
  double keff_diff = 100.;
  double outer_flux_diff = 100.;
  std::size_t outer_iter = 0;
  while (keff_diff > 1.E-5 || outer_flux_diff > 1.E-5) {
    outer_iter++;
    old_outer_flux = flux_;
    fill_fission_source(Qfiss, flux_);

    inner_flux = flux_;
    double inner_flux_diff = 100.; 
    while (inner_flux_diff > 1.E-5) {
      fill_scatter_source(Qscat, inner_flux);

      Q = Qfiss + Qscat;
      //std::cout << "Q = " << Q << "\n";

      sweep(next_inner_flux, Q);

      // Get difference
      inner_flux_diff = xt::amax(xt::abs(next_inner_flux - inner_flux) / next_inner_flux)();
      //std::cout << ">>> Inner it diff = " << inner_flux_diff << "\n";
      //char v;
      //std::cin >> v;

      inner_flux = next_inner_flux;
      next_inner_flux.fill(0.);
    }
    flux_ = inner_flux;
    
    // Get difference
    outer_flux_diff = xt::amax(xt::abs(old_outer_flux - flux_) / flux_)();
    const double old_keff = keff_;
    keff_ = calc_keff(old_outer_flux, flux_, keff_);
    keff_diff = std::abs((old_keff - keff_)/keff_);

    spdlog::info("-------------------------------------");
    spdlog::info("Iteration {:>4d}          keff: {:.5f}", outer_iter, keff_);
    spdlog::info("     keff difference:     {:.5E}", keff_diff);
    spdlog::info("     max flux difference: {:.5E}", outer_flux_diff);
  }
}

void ReflectorSN::sweep(xt::xtensor<double, 2>& flux, const xt::xtensor<double, 2>& Q) {
  const std::size_t NG = xs_.front()->ngroups();

  xt::xtensor<double, 1> incident_angular_flux = xt::zeros<double>({mu_.size()});

  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t n = 0; n < mu_.size(); n++) {
      const double mu = mu_[n];
      double flux_in = 0.;
      double flux_out = 0.;
      double flux_bin = 0.;

      if (mu < 0.) {
        // Track from right to left (negative direction)
        flux_in = 0.;

        for (int ii = static_cast<int>(xs_.size())-1; ii >= 0; ii--) {
          const std::size_t i = static_cast<std::size_t>(ii);
          const double dx = dx_[i];
          const double Etr = xs_[i]->Etr(g);
          const double Qni = Q(g, i);

          // Calculate outgoing flux and bin flux
          flux_out = (2.*dx*Qni + (2.*std::abs(mu) - dx*Etr)*flux_in) / (dx*Etr + 2.*std::abs(mu));
          flux_bin = 0.5*(flux_in + flux_out);

          // Contribute to flux legendre moments
          flux(g,i) += wgt_[n] * flux_bin;

          // Save outgoing flux as an incident flux
          if (i == 0) {
            incident_angular_flux[mu_.size()-1-n] = flux_out;
          }

          flux_in = flux_out;
        }
      } else {
        // Track from left to right (positive direction)
        flux_in = incident_angular_flux[n];

        for (std::size_t i = 0; i < xs_.size(); i++) {
          const double dx = dx_[i];
          const double Etr = xs_[i]->Etr(g);
          const double Qni = Q(g, i);

          // Calculate outgoing flux and bin flux
          flux_out = (2.*dx*Qni + (2.*std::abs(mu) - dx*Etr)*flux_in) / (dx*Etr + 2.*std::abs(mu));
          flux_bin = 0.5*(flux_in + flux_out);

          // Contribute to flux legendre moments
          flux(g,i) += wgt_[n] * flux_bin;

          flux_in = flux_out;
        }
      }
    } // for all mu
    incident_angular_flux.fill(0.);
  } // for all groups
}

double ReflectorSN::calc_keff(const xt::xtensor<double, 2>& old_flux, const xt::xtensor<double, 2>& new_flux, const double keff) const {
  double num = 0.;
  double denom = 0.;
  for (std::size_t i = 0; i < xs_.size(); i++) {
    const auto& mat = xs_[i];
    const double dx = dx_[i];
    for (std::size_t g = 0; g < mat->ngroups(); g++) {
      num += dx * mat->vEf(g) * new_flux(g, i);
      denom += dx * mat->vEf(g) * old_flux(g, i);
    }
  }

  return keff * num / denom;
}

void ReflectorSN::fill_fission_source(xt::xtensor<double, 2>& Qfiss, const xt::xtensor<double, 2>& flux) const {
  Qfiss.fill(0.);

  for (std::size_t i = 0; i < xs_.size(); i++) {
    const auto& mat = xs_[i];

    for (std::size_t g = 0; g < xs_[0]->ngroups(); g++) {
      const double chi_g = mat->chi(g);

      for (std::size_t gg = 0; gg < xs_[0]->ngroups(); gg++) {
        Qfiss(g, i) += chi_g * mat->vEf(gg) * flux(gg, i);
      }
    }
  }

  Qfiss /= 2. * keff_;
}

void ReflectorSN::fill_scatter_source(xt::xtensor<double, 2>& Qscat, const xt::xtensor<double, 2>& flux) const {
  Qscat.fill(0.);

  for (std::size_t i = 0; i < xs_.size(); i++) {
    const auto& mat = xs_[i];

    for (std::size_t g = 0; g < xs_[0]->ngroups(); g++) {
      for (std::size_t gg = 0; gg < xs_[0]->ngroups(); gg++) {
        Qscat(g, i) += 0.5 * mat->Es_tr(gg, g) * flux(gg, i);
      }
    }
  }
}

}