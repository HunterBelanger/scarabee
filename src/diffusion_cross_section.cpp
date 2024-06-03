#include <diffusion_cross_section.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xbuilder.hpp>

#include <cmath>
#include <sstream>

namespace scarabee {

DiffusionCrossSection::DiffusionCrossSection(const xt::xtensor<double, 1>& D,
                                             const xt::xtensor<double, 1>& Ea,
                                             const xt::xtensor<double, 2>& Es,
                                             const xt::xtensor<double, 1>& Ef,
                                             const xt::xtensor<double, 1>& vEf,
                                             const xt::xtensor<double, 1>& chi,
                                             const std::string& name)
    : Es_(Es),
      D_(D),
      Ea_(Ea),
      Ef_(Ef),
      vEf_(vEf),
      chi_(chi),
      name_(name),
      fissile_(false) {
  this->check_xs();
}

DiffusionCrossSection::DiffusionCrossSection(const xt::xtensor<double, 1>& D,
                                             const xt::xtensor<double, 1>& Ea,
                                             const xt::xtensor<double, 2>& Es,
                                             const std::string& name)
    : Es_(Es), D_(D), Ea_(Ea), vEf_(), chi_(), name_(name), fissile_(false) {
  Ef_ = xt::zeros<double>({ngroups()});
  vEf_ = xt::zeros<double>({ngroups()});
  chi_ = xt::zeros<double>({ngroups()});
  this->check_xs();
}

void DiffusionCrossSection::check_xs() {
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

  if (Es_.shape(0) != ngroups() || Es_.shape(1) != ngroups()) {
    auto mssg = "Es is not the same size of Et.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check values for each group
  double chi_sum = 0.;
  for (std::size_t g = 0; g < ngroups(); g++) {
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

    if (std::isnan(D_(gin))) {
      std::stringstream mssg;
      mssg << "D has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (D_(gin) == 0.) {
      std::stringstream mssg;
      mssg << "D value of zero in group " << gin << ".";
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
      if (std::isnan(Es_(gin, gout))) {
        std::stringstream mssg;
        mssg << "Es has NaN value in transfer " << gin << " -> " << gout << ".";
        spdlog::error(mssg.str());
        throw ScarabeeException(mssg.str());
      }
    }
  }
}

}  // namespace scarabee
