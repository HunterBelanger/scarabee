#include <mg_cross_sections.hpp>

MGCrossSections::MGCrossSections() {}

MGCrossSections& MGCrossSections::operator+=(const MGCrossSections& R) {
  // Perform dimensionality checks
  if (ngroups() != R.ngroups()) {
    throw ScarabeeException("Disagreement in number of groups.");
  }

  if (ndelayed_families() != R.ndelayed_families()) {
    throw ScarabeeException("Disagreement in number of delayed families.");
  }

  if (max_legendre_order() != R.max_legendre_order()) {
    throw ScarabeeException("Disagreement in maximum legendre order.");
  }

  // First, calculate/update fission quantities which must be averaged
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    // Get new fissioin xs
    double Ef_g = Ef(g) + R.Ef(g);

    // Calc new prompt data
    double nu_p_g = (nu_prompt(g) * Ef(g) + R.nu_prompt(g) * R.Ef(g)) / Ef_g;
    double chi_p_g = (chi_prompt(g) * nu_prompt(g) * Ef(g) +
                      R.chi_prompt(g) * R.nu_prompt(g) * R.Ef(g)) /
                     (nu_p_g * Ef_g);

    // Assign new prompt data
    nu_(0, g) = nu_p_g;
    chi_(0, g) = chi_p_g;

    // Get new delayed data
    for (std::uint32_t i = 1; i <= ndelayed_families(); i++) {
      // Calc new prompt data
      double nu_i_g =
          (nu_delayed(g, i) * Ef(g) + R.nu_delayed(g, i) * R.Ef(g)) / Ef_g;
      double chi_i_g = (chi_delayed(g, i) * nu_delayed(g, i) * Ef(g) +
                        R.chi_delayed(g, i) * R.nu_delayed(g, i) * R.Ef(g)) /
                       (nu_i_g * Ef_g);

      // Assign new prompt data
      nu_(i, g) = nu_i_g;
      chi_(i, g) = chi_i_g;
    }

    Ef_[g] = Ef_g;
  }

  // Make sure that we are fissile if the other was also fissile
  if (R.fissile()) fissile_ = true;

  // Add all other cross sections which aren't averaged
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    for (std::uint32_t gout = 0; gout < ngroups(); gout++) {
      for (std::uint32_t l = 0; l <= max_legendre_order(); l++) {
        Es_(g, gout, l) += R.Es(g, gout, l);
      }

      Es_tr_(g, gout) += R.Es_tr(g, gout);
    }

    Et_[g] += R.Et(g);
    Etr_[g] += R.Etr(g);
    Ea_[g] += R.Ea(g);
    Ef_[g] += R.Ef(g);
  }

  return *this;
}

MGCrossSections& MGCrossSections::operator*=(double N) {
  // Scale Es_
  for (auto& v : Es_.data_vector()) v *= N;

  // Scale Es_tr_
  for (auto& v : Es_tr_.data_vector()) v *= N;

  // Scale Et_
  for (auto& v : Et_) v *= N;

  // Scale Etr_
  for (auto& v : Etr_) v *= N;

  // Scale Ea_
  for (auto& v : Ea_) v *= N;

  // Scale Ef_
  for (auto& v : Ef_) v *= N;

  return *this;
}

MGCrossSections MGCrossSections::operator+(const MGCrossSections& R) const {
  MGCrossSections out = *this;
  out += R;
  return out;
}

MGCrossSections MGCrossSections::operator*(double N) const {
  MGCrossSections out = *this;
  out *= N;
  return out;
}

void MGCrossSections::generate_transport_corrected_data() {
  // Now make a modified scattering matrix
  Es_tr_.reshape({ngroups(), ngroups()});
  Etr_.resize(ngroups());

  // Now go through all groups
  for (std::uint32_t g = 0; g < ngroups(); g++) {
    // Go through all out groups, and make transport corrected scatter matrix
    for (std::uint32_t gout = 0; gout < ngroups(); gout++) {
      Es_tr_(g, gout) = Es(g, gout, 0);

      if (g == gout && max_legendre_order() > 0) {
        Es_tr_(g, g) -= Es(g, g, 1);
      }
    }

    // Make transport cross section
    if (max_legendre_order() > 0)
      Etr_[g] = Et(g) - Es(g, g, 1);
    else
      Etr_[g] = Et(g);
  }
}
