#include <cross_section.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xbuilder.hpp>

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
