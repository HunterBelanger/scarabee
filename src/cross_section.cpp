#include <cross_section.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xbuilder.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace scarabee {

CrossSection::CrossSection(const xt::xtensor<double, 1>& Etr,
                           const xt::xtensor<double, 1>& Ea,
                           const xt::xtensor<double, 2>& Es_tr,
                           const xt::xtensor<double, 1>& Ef,
                           const xt::xtensor<double, 1>& vEf,
                           const xt::xtensor<double, 1>& chi,
                           const std::string& name)
    : Es_tr_(Es_tr),
      Es1_(),
      Etr_(Etr),
      Ea_(Ea),
      Ef_(Ef),
      vEf_(vEf),
      chi_(chi),
      name_(name),
      fissile_(false) {
  this->check_xs();
}

CrossSection::CrossSection(
    const xt::xtensor<double, 1>& Et, const xt::xtensor<double, 1>& Ea,
    const xt::xtensor<double, 2>& Es, const xt::xtensor<double, 2>& Es1,
    const xt::xtensor<double, 1>& Ef, const xt::xtensor<double, 1>& vEf,
    const xt::xtensor<double, 1>& chi, const std::string& name)
    : Es_tr_(Es),
      Es1_(Es1),
      Etr_(Et),
      Ea_(Ea),
      Ef_(Ef),
      vEf_(vEf),
      chi_(chi),
      name_(name),
      fissile_(false) {
  // Check xs first, to make sure the size of everything is okay
  this->check_xs();

  // Apply the transport correction
  for (std::size_t g = 0; g < ngroups(); g++) {
    Etr_(g) -= Es1_(g, g);
    Es_tr_(g, g) -= Es1_(g, g);
  }
}

CrossSection::CrossSection(const xt::xtensor<double, 1>& Etr,
                           const xt::xtensor<double, 1>& Ea,
                           const xt::xtensor<double, 2>& Es_tr,
                           const std::string& name)
    : Es_tr_(Es_tr),
      Es1_(),
      Etr_(Etr),
      Ea_(Ea),
      vEf_(),
      chi_(),
      name_(name),
      fissile_(false) {
  Ef_ = xt::zeros<double>({ngroups()});
  vEf_ = xt::zeros<double>({ngroups()});
  chi_ = xt::zeros<double>({ngroups()});
  this->check_xs();
}

CrossSection::CrossSection(const xt::xtensor<double, 1>& Et,
                           const xt::xtensor<double, 1>& Ea,
                           const xt::xtensor<double, 2>& Es,
                           const xt::xtensor<double, 2>& Es1,
                           const std::string& name)
    : Es_tr_(Es),
      Es1_(Es1),
      Etr_(Et),
      Ea_(Ea),
      vEf_(),
      chi_(),
      name_(name),
      fissile_(false) {
  Ef_ = xt::zeros<double>({ngroups()});
  vEf_ = xt::zeros<double>({ngroups()});
  chi_ = xt::zeros<double>({ngroups()});
  this->check_xs();

  // Apply the transport correction
  for (std::size_t g = 0; g < ngroups(); g++) {
    Etr_(g) -= Es1_(g, g);
    Es_tr_(g, g) -= Es1_(g, g);
  }
}

std::shared_ptr<CrossSection> CrossSection::condense(const std::vector<std::pair<std::size_t, std::size_t>>& groups, const std::vector<double>& flux) const {
  const std::size_t NG = ngroups();

  // First, we need to check the condensation scheme
  if (groups.size() == 0) {
    auto mssg = "Empty energy condensation scheme provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (groups.front().first != 0) {
    auto mssg = "The energy condensation scheme does not start with 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (groups.back().second != NG-1) {
    std::stringstream mssg;
    mssg << "The energy condensation scheme does not end with " << NG-1 << ".";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  for (std::size_t i = 0; i < NG-1; i++) {
    if (groups[i].second + 1 != groups[i+1].first) {
      std::stringstream mssg;
      mssg << "Condensed groups " << i << " and " << i+1 << " are not continuous.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  if (flux.size() != NG) {
    auto mssg = "The number of provided flux values diagrees with the number of energy groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  
  // Everything checks out, we can start the condensation.

  const std::size_t NGOUT = groups.size();
  xt::xtensor<double, 1> Et  = xt::zeros<double>({NGOUT});
  xt::xtensor<double, 1> Ea  = xt::zeros<double>({NGOUT});
  xt::xtensor<double, 2> Es  = xt::zeros<double>({NGOUT, NGOUT});
  xt::xtensor<double, 1> Ef  = xt::zeros<double>({NGOUT});
  xt::xtensor<double, 1> vEf = xt::zeros<double>({NGOUT});
  xt::xtensor<double, 1> chi = xt::zeros<double>({NGOUT});
  xt::xtensor<double, 2> Es1;

  const bool has_P1 = this->anisotropic();
  if (has_P1) Es1 = xt::zeros<double>({NGOUT, NGOUT});

  for (std::size_t G = 0; G < NGOUT; G++) { // Incoming macro groups
    const std::size_t g_min = groups[G].first;
    const std::size_t g_max = groups[G].second;

    // First, we get the sum of all flux values in the macro group
    const double flux_G = std::accumulate(flux.begin()+g_min, flux.begin()+g_max, 0.);
    const double invs_flux_G = 1. / flux_G;

    // First we do all of the 1D cross sections
    for (std::size_t g = g_min; g <= g_max; g++) {
      Et(G)  += invs_flux_G * flux[g] * this->Et(g);
      Ea(G)  += invs_flux_G * flux[g] * this->Ea(g);
      Ea(G)  += invs_flux_G * flux[g] * this->Ea(g);
      Ef(G)  += invs_flux_G * flux[g] * this->Ef(g);
      vEf(G) += invs_flux_G * flux[g] * this->vEf(g);
      chi(G) += this->chi(g); // chi doesn't need to be weighted
    }

    // Next, we do all 2D cross sections
    for (std::size_t GG = 0; GG < NGOUT; GG++) { // Outgoing macro groups
      const std::size_t gg_min = groups[GG].first;
      const std::size_t gg_max = groups[GG].second;
      for (std::size_t g = g_min; g <= g_max; g++) { // Incoming micro groups
        for (std::size_t gg = gg_min; gg <= gg_max; gg++) { // Outgoing micro groups
          Es(G, GG) += invs_flux_G * this->Es(g, gg);

          if (has_P1) Es1(G, GG) += invs_flux_G * this->Es1(g, gg);
        }
      }
    }
  }

  if (has_P1) {
    return std::make_shared<CrossSection>(Et, Ea, Es, Es1, Ef, vEf, chi);
  }

  // If we don't have a P1 matrix, then Et actually contains Etr and Es
  // contains Es_tr. We can therefore use that constructor.
  return std::make_shared<CrossSection>(Et, Ea, Es, Ef, vEf, chi);
}

CrossSection& CrossSection::operator+=(const CrossSection& R) {
  // Perform dimensionality checks
  if (ngroups() != R.ngroups()) {
    auto mssg = "Disagreement in number of groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, calculate/update fission spectrum
  double L_vEf_sum = 0.;
  double R_vEf_sum = 0.;
  for (std::size_t g = 0; g < ngroups(); g++) {
    L_vEf_sum += vEf(g);
    R_vEf_sum += R.vEf(g);
  }

  double chi_sum = 0.;
  if (L_vEf_sum + R_vEf_sum > 0.) {
    for (std::size_t g = 0; g < ngroups(); g++) {
      chi_(g) =
          (chi(g) * L_vEf_sum + R.chi(g) * R_vEf_sum) / (L_vEf_sum + R_vEf_sum);
      chi_sum += chi_(g);
    }

    // Renormalize chi
    if (chi_sum > 0.) chi_ /= chi_sum;
  } else {
    chi_.fill(0.);
  }

  // Make sure that we are fissile if the other was also fissile
  if (R.fissile()) fissile_ = true;

  // Add all other cross sections which aren't averaged
  for (std::size_t g = 0; g < ngroups(); g++) {
    for (std::size_t gout = 0; gout < ngroups(); gout++) {
      Es_tr_(g, gout) += R.Es_tr(g, gout);
    }

    Etr_(g) += R.Etr(g);
    Ea_(g) += R.Ea(g);
    Ef_(g) += R.Ef(g);
    vEf_(g) += R.vEf(g);
  }

  this->check_xs();

  return *this;
}

CrossSection& CrossSection::operator*=(double N) {
  // Scale Es_tr_
  Es_tr_ *= N;

  // Scale Etr_
  Etr_ *= N;

  // Scale Ea_
  Ea_ *= N;

  // Scale Ef_
  Ef_ *= N;

  // Scale vEf_
  vEf_ *= N;

  this->check_xs();

  return *this;
}

CrossSection CrossSection::operator+(const CrossSection& R) const {
  CrossSection out = *this;
  out += R;
  return out;
}

CrossSection CrossSection::operator*(double N) const {
  CrossSection out = *this;
  out *= N;
  return out;
}

void CrossSection::check_xs() {
  // Make correct number of groups
  if (ngroups() == 0) {
    auto mssg = "Must have at least 1 energy group.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Ea_.size() != ngroups()) {
    auto mssg = "Ea is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Ef_.size() != ngroups()) {
    auto mssg = "Ef is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (vEf_.size() != ngroups()) {
    auto mssg = "vEf is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (chi_.size() != ngroups()) {
    auto mssg = "chi is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Es_tr_.shape(0) != ngroups() || Es_tr_.shape(1) != ngroups()) {
    auto mssg = "Es is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Es1_.size() > 0) {
    if (Es1_.shape(0) != ngroups() || Es1_.shape(1) != ngroups()) {
      auto mssg = "Es1 is not the same size of Et.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Check values for each group
  double chi_sum = 0.;
  for (std::size_t g = 0; g < ngroups(); g++) {
    if (std::abs((Etr(g) - Ea(g) - Es_tr(g)) / Etr(g)) > 1.E-3) {
      std::stringstream mssg;
      mssg << "Ea + Es_tr != Etr in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (vEf(g) < 0.) {
      std::stringstream mssg;
      mssg << "vEf is negative in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (vEf(g) > 0.) fissile_ = true;

    if (chi(g) < 0.) {
      std::stringstream mssg;
      mssg << "chi is negative in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
    chi_sum += chi(g);
  }

  if (fissile_ && chi_sum <= 0.) {
    std::stringstream mssg;
    mssg << "Material is fissile, but all chi in all groups is zero.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Make sure no NaN values
  for (std::size_t gin = 0; gin < ngroups(); gin++) {
    if (std::isnan(Ea_(gin))) {
      std::stringstream mssg;
      mssg << "Ea has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (std::isnan(Ef_(gin))) {
      std::stringstream mssg;
      mssg << "Ef has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (std::isnan(Etr_(gin))) {
      std::stringstream mssg;
      mssg << "Etr has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (Etr_(gin) == 0.) {
      std::stringstream mssg;
      mssg << "Etr value of zero in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (std::isnan(vEf_(gin))) {
      std::stringstream mssg;
      mssg << "vEf has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (std::isnan(chi_(gin))) {
      std::stringstream mssg;
      mssg << "chi has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    for (std::size_t gout = 0; gout < ngroups(); gout++) {
      if (std::isnan(Es_tr_(gin, gout))) {
        std::stringstream mssg;
        mssg << "Es_tr has NaN value in transfer " << gin << " -> " << gout
             << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }

      if (Es1_.size() > 0 && std::isnan(Es1_(gin, gout))) {
        std::stringstream mssg;
        mssg << "Es1 has NaN value in transfer " << gin << " -> " << gout
             << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }
    }
  }
}

}  // namespace scarabee
