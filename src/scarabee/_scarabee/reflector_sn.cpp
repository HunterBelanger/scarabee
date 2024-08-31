#include <reflector_sn.hpp>
#include <utils/math.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xio.hpp>

#include <sstream>

namespace scarabee {

ReflectorSN::ReflectorSN(const std::vector<std::shared_ptr<CrossSection>>& xs,
                         const xt::xtensor<double, 1>& dx)
    : xs_(xs), dx_(dx) {
  if (xs_.size() != dx_.size()) {
    auto mssg = "Number of cross sections and regions do not agree.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (xs_.size() == 0) {
    auto mssg = "Number of regions must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < dx_.size(); i++) {
    if (dx_[i] <= 0.) {
      std::stringstream mssg;
      mssg << "Region " << i << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  ngroups_ = xs_.front()->ngroups();
  for (std::size_t i = 0; i < xs_.size(); i++) {
    if (xs_[i]->ngroups() != ngroups_) {
      auto mssg = "Not all regions have the same number of groups.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Must allocate with zeros in case someone calls the flux method
  flux_ = xt::zeros<double>({ngroups_, xs_.size()});
}

void ReflectorSN::set_flux_tolerance(double ftol) {
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

void ReflectorSN::set_keff_tolerance(double ktol) {
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
  while (keff_diff > keff_tol_ || outer_flux_diff > flux_tol_) {
    outer_iter++;
    old_outer_flux = flux_;
    fill_fission_source(Qfiss, flux_);

    inner_flux = flux_;
    double inner_flux_diff = 100.;
    while (inner_flux_diff > flux_tol_) {
      fill_scatter_source(Qscat, inner_flux);

      Q = Qfiss + Qscat;

      sweep(next_inner_flux, Q);

      // Get difference
      inner_flux_diff =
          xt::amax(xt::abs(next_inner_flux - inner_flux) / next_inner_flux)();

      inner_flux = next_inner_flux;
      next_inner_flux.fill(0.);
    }
    flux_ = inner_flux;

    // Get difference
    outer_flux_diff = xt::amax(xt::abs(old_outer_flux - flux_) / flux_)();
    const double old_keff = keff_;
    keff_ = calc_keff(old_outer_flux, flux_, keff_);
    keff_diff = std::abs((old_keff - keff_) / keff_);

    spdlog::info("-------------------------------------");
    spdlog::info("Iteration {:>4d}          keff: {:.5f}", outer_iter, keff_);
    spdlog::info("     keff difference:     {:.5E}", keff_diff);
    spdlog::info("     max flux difference: {:.5E}", outer_flux_diff);
  }

  solved_ = true;
}

void ReflectorSN::sweep(xt::xtensor<double, 2>& flux,
                        const xt::xtensor<double, 2>& Q) {
  const std::size_t NG = xs_.front()->ngroups();

  xt::xtensor<double, 1> incident_angular_flux =
      xt::zeros<double>({mu_.size()});

  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t n = 0; n < mu_.size(); n++) {
      const double mu = mu_[n];
      double flux_in = 0.;
      double flux_out = 0.;
      double flux_bin = 0.;

      if (mu < 0.) {
        // Track from right to left (negative direction)
        flux_in = 0.;

        for (int ii = static_cast<int>(xs_.size()) - 1; ii >= 0; ii--) {
          const std::size_t i = static_cast<std::size_t>(ii);
          const double dx = dx_[i];
          const double Etr = xs_[i]->Etr(g);
          const double Qni = Q(g, i);

          // Calculate outgoing flux and bin flux
          flux_out =
              (2. * dx * Qni + (2. * std::abs(mu) - dx * Etr) * flux_in) /
              (dx * Etr + 2. * std::abs(mu));
          flux_bin = 0.5 * (flux_in + flux_out);

          // Contribute to flux legendre moments
          flux(g, i) += wgt_[n] * flux_bin;

          // Save outgoing flux as an incident flux
          if (i == 0) {
            incident_angular_flux[mu_.size() - 1 - n] = flux_out;
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
          flux_out =
              (2. * dx * Qni + (2. * std::abs(mu) - dx * Etr) * flux_in) /
              (dx * Etr + 2. * std::abs(mu));
          flux_bin = 0.5 * (flux_in + flux_out);

          // Contribute to flux legendre moments
          flux(g, i) += wgt_[n] * flux_bin;

          flux_in = flux_out;
        }
      }
    }  // for all mu
    incident_angular_flux.fill(0.);
  }  // for all groups
}

double ReflectorSN::calc_keff(const xt::xtensor<double, 2>& old_flux,
                              const xt::xtensor<double, 2>& new_flux,
                              const double keff) const {
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

void ReflectorSN::fill_fission_source(
    xt::xtensor<double, 2>& Qfiss, const xt::xtensor<double, 2>& flux) const {
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

void ReflectorSN::fill_scatter_source(
    xt::xtensor<double, 2>& Qscat, const xt::xtensor<double, 2>& flux) const {
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

double ReflectorSN::flux(std::size_t i, std::size_t g) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "Region index i =" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Energy group index g =" << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return flux_(g, i);
}

const std::shared_ptr<CrossSection> ReflectorSN::xs(std::size_t i) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "Region index i =" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return xs_[i];
}

double ReflectorSN::volume(std::size_t i) const {
  if (i >= this->size()) {
    std::stringstream mssg;
    mssg << "Region index i =" << i << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return dx_[i];
}

std::shared_ptr<CrossSection> ReflectorSN::homogenize(
    const std::vector<std::size_t>& regions) const {
  // We can only perform a homogenization if we have a flux spectrum
  if (solved() == false) {
    auto mssg =
        "Cannot perform homogenization when problem has not been solved.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check all regions are valid
  if (regions.size() > this->nregions()) {
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

  // We now begin homogenization
  const std::size_t NR = regions.size();
  const std::size_t NG = this->ngroups();

  bool has_P1 = false;
  for (const auto m : regions) {
    if (this->xs(m)->anisotropic()) {
      has_P1 = true;
      break;
    }
  }

  xt::xtensor<double, 1> Et = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 2> Es = xt::zeros<double>({NG, NG});
  xt::xtensor<double, 1> Ef = xt::zeros<double>({NG});
  xt::xtensor<double, 1> vEf = xt::zeros<double>({NG});
  xt::xtensor<double, 1> chi = xt::zeros<double>({NG});
  xt::xtensor<double, 2> Es1;
  if (has_P1) Es1 = xt::zeros<double>({NG, NG});

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
      sum_fluxV += this->flux(i, g) * dx_[i];
    }
    const double invs_sum_fluxV = 1. / sum_fluxV;

    j = 0;
    for (const auto i : regions) {
      const auto& mat = this->xs(i);
      const double V = this->volume(i);
      const double flx = flux(i, g);
      const double coeff = invs_sum_fluxV * flx * V;
      Ea(g) += coeff * mat->Ea(g);
      Ef(g) += coeff * mat->Ef(g);
      vEf(g) += coeff * mat->vEf(g);

      chi(g) += invs_sum_fiss_prod * fiss_prod[j] * mat->chi(g);

      for (std::size_t gg = 0; gg < NG; gg++) {
        Es(g, gg) += coeff * mat->Es(g, gg);

        if (has_P1) Es1(g, gg) += coeff * mat->Es1(g, gg);
      }

      j++;
    }

    // Reconstruct total xs from absorption and scattering
    Et(g) = Ea(g) + xt::sum(xt::view(Es, g, xt::all()))();
  }

  if (has_P1) {
    return std::make_shared<CrossSection>(Et, Ea, Es, Es1, Ef, vEf, chi);
  }

  // If we don't have a P1 matrix, then Et actually contains Etr and Es
  // contains Es_tr. We can therefore use that constructor.
  return std::make_shared<CrossSection>(Et, Ea, Es, Ef, vEf, chi);
}

xt::xtensor<double, 1> ReflectorSN::homogenize_flux_spectrum(
    const std::vector<std::size_t>& regions) const {
  // We can only perform a homogenization if we have a flux spectrum
  if (solved() == false) {
    auto mssg =
        "Cannot perform spectrum homogenization when problem has not been "
        "solved.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check all regions are valid
  if (regions.size() > this->nregions()) {
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

const std::array<double, 32> ReflectorSN::mu_{
    -9.97263861849481563545e-01, -9.85611511545268335400e-01,
    -9.64762255587506430774e-01, -9.34906075937739689171e-01,
    -8.96321155766052123965e-01, -8.49367613732569970134e-01,
    -7.94483795967942406963e-01, -7.32182118740289680387e-01,
    -6.63044266930215200975e-01, -5.87715757240762329041e-01,
    -5.06899908932229390024e-01, -4.21351276130635345364e-01,
    -3.31868602282127649780e-01, -2.39287362252137074545e-01,
    -1.44471961582796493485e-01, -4.83076656877383162348e-02,
    4.83076656877383162348e-02,  1.44471961582796493485e-01,
    2.39287362252137074545e-01,  3.31868602282127649780e-01,
    4.21351276130635345364e-01,  5.06899908932229390024e-01,
    5.87715757240762329041e-01,  6.63044266930215200975e-01,
    7.32182118740289680387e-01,  7.94483795967942406963e-01,
    8.49367613732569970134e-01,  8.96321155766052123965e-01,
    9.34906075937739689171e-01,  9.64762255587506430774e-01,
    9.85611511545268335400e-01,  9.97263861849481563545e-01};

const std::array<double, 32> ReflectorSN::wgt_{
    7.01861000947009660041e-03, 1.62743947309056706052e-02,
    2.53920653092620594558e-02, 3.42738629130214331027e-02,
    4.28358980222266806569e-02, 5.09980592623761761962e-02,
    5.86840934785355471453e-02, 6.58222227763618468377e-02,
    7.23457941088485062254e-02, 7.81938957870703064717e-02,
    8.33119242269467552222e-02, 8.76520930044038111428e-02,
    9.11738786957638847129e-02, 9.38443990808045656392e-02,
    9.56387200792748594191e-02, 9.65400885147278005668e-02,
    9.65400885147278005668e-02, 9.56387200792748594191e-02,
    9.38443990808045656392e-02, 9.11738786957638847129e-02,
    8.76520930044038111428e-02, 8.33119242269467552222e-02,
    7.81938957870703064717e-02, 7.23457941088485062254e-02,
    6.58222227763618468377e-02, 5.86840934785355471453e-02,
    5.09980592623761761962e-02, 4.28358980222266806569e-02,
    3.42738629130214331027e-02, 2.53920653092620594558e-02,
    1.62743947309056706052e-02, 7.01861000947009660041e-03};
}  // namespace scarabee