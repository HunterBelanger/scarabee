#include <transport_xs.hpp>

TransportXS::TransportXS() {}

TransportXS& TransportXS::operator+=(const TransportXS& R) {
  // Perform dimensionality checks
  if (ngroups() != R.ngroups()) {
    throw ScarabeeException("Disagreement in number of groups.");
  }

  // First, calculate/update fission quantities which must be averaged
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    // Get new fissioin xs
    double Ef_g = Ef(g) + R.Ef(g);

    // Calc new fission data
    double nu_p_g = (nu(g) * Ef(g) + R.nu(g) * R.Ef(g)) / Ef_g;
    double chi_p_g = (chi(g) * nu(g) * Ef(g) + R.chi(g) * R.nu(g) * R.Ef(g)) /
                     (nu_p_g * Ef_g);

    // Assign new prompt data
    nu_(g) = nu_p_g;
    chi_(g) = chi_p_g;
    Ef_(g) = Ef_g;
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
    Ef_(g) += R.Ef(g);
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

  // Scale Ef_
  Ef_ *= N;

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
