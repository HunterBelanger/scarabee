#include <transport_xs.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xbuilder.hpp>

#include <cmath>
#include <sstream>

namespace scarabee {

TransportXS::TransportXS(const xt::xtensor<double, 1>& Et,
                         const xt::xtensor<double, 1>& Ea,
                         const xt::xtensor<double, 2>& Es,
                         const xt::xtensor<double, 1>& Ef,
                         const xt::xtensor<double, 1>& vEf,
                         const xt::xtensor<double, 1>& chi,
                         const std::string& name)
    : Es_(Es),
      Et_(Et),
      Ea_(Ea),
      Ef_(Ef),
      vEf_(vEf),
      chi_(chi),
      name_(name),
      fissile_(false) {
  this->check_xs();
}

TransportXS::TransportXS(const xt::xtensor<double, 1>& Et,
                         const xt::xtensor<double, 1>& Ea,
                         const xt::xtensor<double, 2>& Es,
                         const std::string& name)
    : Es_(Es), Et_(Et), Ea_(Ea), vEf_(), chi_(), name_(name), fissile_(false) {
  Ef_ = xt::zeros<double>({ngroups()});
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
  double L_vEf_sum = 0.;
  double R_vEf_sum = 0.;
  for (std::size_t g = 0; g < ngroups(); g++) {
    L_vEf_sum += vEf(g);
    R_vEf_sum += R.vEf(g);
  }

  double chi_sum = 0.; 
  if (L_vEf_sum+R_vEf_sum > 0.) {
    for (std::size_t g = 0; g < ngroups(); g++) {
      chi_(g) = (chi(g)*L_vEf_sum + R.chi(g)*R_vEf_sum) / (L_vEf_sum + R_vEf_sum);
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
      Es_(g, gout) += R.Es(g, gout);
    }

    Et_(g) += R.Et(g);
    Ea_(g) += R.Ea(g);
    Ef_(g) += R.Ef(g);
    vEf_(g) += R.vEf(g);
  }

  this->check_xs();

  return *this;
}

TransportXS& TransportXS::operator*=(double N) {
  // Scale Es_
  Es_ *= N;

  // Scale Et_
  Et_ *= N;

  // Scale Ea_
  Ea_ *= N;

  // Scale Ef_
  Ef_ *= N;

  // Scale vEf_
  vEf_ *= N;

  this->check_xs();

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
    if (std::abs((Et(g) - Ea(g) - Es(g)) / Et(g)) > 1.E-3) {
      std::stringstream mssg;
      mssg << "Ea + Es != Et in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    /*
    if (Ea(g) < 0.) {
      std::stringstream mssg;
      mssg << "Ea is negative in group " << g << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
    */

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

    if (std::isnan(Et_(gin))) {
      std::stringstream mssg;
      mssg << "Et has NaN value in group " << gin << ".";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (Et_(gin) == 0.) {
      std::stringstream mssg;
      mssg << "Et value of zero in group " << gin << ".";
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
