#include <transport_xs.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xbuilder.hpp>

#include <sstream>

namespace scarabee {

TransportXS::TransportXS(const xt::xtensor<double, 1>& Et,
                         const xt::xtensor<double, 1>& Ea,
                         const xt::xtensor<double, 2>& Es,
                         const xt::xtensor<double, 1>& vEf,
                         const xt::xtensor<double, 1>& chi)
    : Es_(Es), Et_(Et), Ea_(Ea), vEf_(vEf), chi_(chi), fissile_(false) {
  this->check_xs();
}

TransportXS::TransportXS(const xt::xtensor<double, 1>& Et,
                         const xt::xtensor<double, 1>& Ea,
                         const xt::xtensor<double, 2>& Es)
    : Es_(Es), Et_(Et), Ea_(Ea), vEf_(), chi_(), fissile_(false) {
  vEf_ = xt::zeros<double>({ngroups()});
  chi_ = xt::zeros<double>({ngroups()});
  this->check_xs();
}

TransportXS& TransportXS::operator+=(const TransportXS& R) {
  // Perform dimensionality checks
  if (ngroups() != R.ngroups()) {
    auto mssg = "Disagreement in number of groups.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, calculate/update fission spectrum
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    // Calc new fission data
    chi_(g) = (chi(g) * vEf(g) + R.chi(g) * R.vEf(g)) / (vEf(g) + R.vEf(g));
  }

  // Make sure that we are fissile if the other was also fissile
  if (R.fissile()) fissile_ = true;

  // Add all other cross sections which aren't averaged
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::uint32_t gout = 0; gout < ngroups(); gout++) {
      Es_(g, gout) += R.Es(g, gout);
    }

    Et_(g) += R.Et(g);
    Ea_(g) += R.Ea(g);
    vEf_(g) += R.vEf(g);
  }

  return *this;
}

TransportXS& TransportXS::operator*=(double N) {
  // Scale Es_
  Es_ *= N;

  // Scale Et_
  Et_ *= N;

  // Scale Ea_
  Ea_ *= N;

  // Scale vEf_
  vEf_ *= N;

  return *this;
}

TransportXS TransportXS::operator+(const TransportXS& R) const {
  TransportXS out = *this;
  out += R;
  return out;
}

TransportXS TransportXS::operator*(double N) const {
  TransportXS out = *this;
  out *= N;
  return out;
}

void TransportXS::check_xs() {
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
    if (std::abs((Et(g) - Ea(g) - Es(g)) / Et(g)) > 1.E-3) {
      std::stringstream mssg;
      mssg << "Ea + Es != Et in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (Ea(g) < 0.) {
      std::stringstream mssg;
      mssg << "Ea is negative in group " << g << ".";
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
}

}  // namespace scarabee
