#include <data/material.hpp>
#include <data/nd_library.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <sstream>

namespace scarabee {

void MaterialComposition::add_nuclide(const std::string& name, double frac) {
  if (frac <= 0.) {
    std::stringstream mssg;
    mssg << "Nuclide \"" << name << "\" given negative or zero fraction.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  nuclides.emplace_back();
  nuclides.back().name = name;
  nuclides.back().fraction = frac;
}

void MaterialComposition::add_nuclide(const Nuclide& nuc) {
  if (nuc.fraction <= 0.) {
    std::stringstream mssg;
    mssg << "Nuclide \"" << nuc.name << "\" has negative or zero fraction.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  nuclides.push_back(nuc);
}

Material::Material(const MaterialComposition& comp, double temp,
                   std::shared_ptr<NDLibrary> ndl)
    : composition_(comp),
      temperature_(temp),
      atoms_per_bcm_(-1.),
      grams_per_cm3_(-1.),
      potential_xs_(0.),
      fissile_(false),
      resonant_(false) {
  // Make sure quantities are positive/valid
  if (temp <= 0.) {
    auto mssg = "Material temperature must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (composition_.nuclides.empty()) {
    auto mssg = "Material must have a composition of at least one component.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ndl == nullptr) {
    auto mssg = "Material must be given NDLibrary.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, we get our density, assuming that it can be computed from the sum
  // of the fractions in the composition
  double frac_sum = 0.;
  for (const auto& c : composition_.nuclides) frac_sum += c.fraction;

  if (composition_.fractions == Fraction::Atoms) {
    atoms_per_bcm_ = frac_sum;
  } else {
    grams_per_cm3_ = frac_sum;
  }

  this->normalize_fractions();

  // Now that fractions have been normalized, we can get the average molar mass
  average_molar_mass_ = this->calc_avg_molar_mass(*ndl);

  // Convert to Atoms fractions if necessary
  if (composition_.fractions == Fraction::Weight) {
    for (auto& c : composition_.nuclides) {
      const auto& nuc = ndl->get_nuclide(c.name);
      c.fraction = c.fraction * average_molar_mass_ / (nuc.awr * N_MASS_AMU);
    }
  }

  // Get the missing density
  if (atoms_per_bcm_ < 0.) {
    atoms_per_bcm_ = (grams_per_cm3_ * N_AVAGADRO) / average_molar_mass_;
  } else {
    grams_per_cm3_ = average_molar_mass_ * atoms_per_bcm_ / N_AVAGADRO;
  }

  // Check fissile and resonant, also get potential_xs
  for (const auto& c : composition_.nuclides) {
    const auto& nuc = ndl->get_nuclide(c.name);
    potential_xs_ += atoms_per_bcm_ * c.fraction * nuc.potential_xs;

    if (nuc.fissile) fissile_ = true;

    if (nuc.resonant) resonant_ = true;
  }
}

Material::Material(const MaterialComposition& comp, double temp, double density,
                   DensityUnits du, std::shared_ptr<NDLibrary> ndl)
    : composition_(comp),
      temperature_(temp),
      atoms_per_bcm_(-1.),
      grams_per_cm3_(-1.),
      potential_xs_(0.),
      fissile_(false),
      resonant_(false) {
  // Make sure quantities are positive/valid
  if (temp <= 0.) {
    auto mssg = "Material temperature must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (density < 0. && du != DensityUnits::sum) {
    auto mssg = "Material density must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (composition_.nuclides.empty()) {
    auto mssg = "Material must have a composition of at least one component.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ndl == nullptr) {
    auto mssg = "Material must be given NDLibrary.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, we get our provided density
  if (du == DensityUnits::sum) {
    double frac_sum = 0.;
    for (const auto& c : composition_.nuclides) frac_sum += c.fraction;

    if (composition_.fractions == Fraction::Atoms) {
      atoms_per_bcm_ = frac_sum;
    } else {
      grams_per_cm3_ = frac_sum;
    }
  } else if (du == DensityUnits::a_bcm) {
    atoms_per_bcm_ = density;
  } else {
    grams_per_cm3_ = density;
  }
  this->normalize_fractions();

  // Now that fractions have been normalized, we can get the average molar mass
  average_molar_mass_ = this->calc_avg_molar_mass(*ndl);

  // Convert to Atoms fractions if necessary
  if (composition_.fractions == Fraction::Weight) {
    for (auto& c : composition_.nuclides) {
      const auto& nuc = ndl->get_nuclide(c.name);
      c.fraction = c.fraction * average_molar_mass_ / (nuc.awr * N_MASS_AMU);
    }
  }

  // Get the missing density
  if (atoms_per_bcm_ < 0.) {
    atoms_per_bcm_ = (grams_per_cm3_ * N_AVAGADRO) / average_molar_mass_;
  } else {
    grams_per_cm3_ = average_molar_mass_ * atoms_per_bcm_ / N_AVAGADRO;
  }

  // Check fissile and resonant, also get potential_xs
  for (const auto& c : composition_.nuclides) {
    const auto& nuc = ndl->get_nuclide(c.name);
    potential_xs_ += atoms_per_bcm_ * c.fraction * nuc.potential_xs;

    if (nuc.fissile) fissile_ = true;

    if (nuc.resonant) resonant_ = true;
  }
}

double Material::atom_density(const std::string& name) const {
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    if (composition_.nuclides[i].name == name) {
      return this->atoms_per_bcm() * composition_.nuclides[i].fraction;
    }
  }

  std::stringstream mssg;
  mssg << "Could not find nuclide with name \"" << name << "\".";
  spdlog::error(mssg.str());
  throw ScarabeeException(mssg.str());

  // NEVER GETS HERE
  return 0.;
}

bool Material::has_component(const std::string& name) const {
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    if (composition_.nuclides[i].name == name) {
      return true;
    }
  }

  return false;
}

double Material::calc_avg_molar_mass(const NDLibrary& ndl) const {
  double avg_mm = 0.;

  for (const auto& comp : composition_.nuclides) {
    const auto& nuc = ndl.get_nuclide(comp.name);
    if (composition_.fractions == Fraction::Atoms) {
      avg_mm += comp.fraction * nuc.awr * N_MASS_AMU;
    } else {
      avg_mm += comp.fraction / (nuc.awr * N_MASS_AMU);
    }
  }

  if (composition_.fractions == Fraction::Weight) avg_mm = 1. / avg_mm;

  return avg_mm;
}

void Material::normalize_fractions() {
  double frac_sum = 0.;

  for (const auto& comp : composition_.nuclides) {
    frac_sum += comp.fraction;
  }

  for (auto& comp : composition_.nuclides) {
    comp.fraction /= frac_sum;
  }
}

std::shared_ptr<CrossSection> Material::carlvik_xs(
    double C, double Ee, std::shared_ptr<NDLibrary> ndl) const {
  // This implementation is based on the methods outlined by Koike and Gibson
  // in his PhD thesis [1,2]. We start by computing the coefficients for the
  // two-term rational approximation, modified according to the Dancoff
  // correction factor, C.
  const double a1 = 0.5 * (C + 5. - std::sqrt(C * C + 34. * C + 1.));
  const double a2 = 0.5 * (C + 5. + std::sqrt(C * C + 34. * C + 1.));
  const double b1 = (a2 - (1. - C)) / (a2 - a1);
  const double b2 = 1. - b1;

  return this->two_term_xs(a1, a2, b1, b2, Ee, ndl);
}

std::shared_ptr<CrossSection> Material::roman_xs(
    double C, double Ee, std::shared_ptr<NDLibrary> ndl) const {
  // This implementation is based on the methods outlined by Koike and Gibson
  // in his PhD thesis [1,2]. We start by computing the coefficients for the
  // two-term rational approximation, modified according to the Dancoff
  // correction factor, C.
  constexpr double a1_base = 1.4;
  constexpr double a2_base = 5.4;
  constexpr double b1_base = 1.1;
  constexpr double b2_base = 0.1;

  if (C > 0.) {
    const double A = (1. - C) / C;
    const double g = A + (b1_base * a1_base) + (b2_base * a2_base);
    const double t = A * (a1_base + a2_base) + a1_base * a2_base;
    const double a1 =
        (t - std::sqrt(t * t - (4. * g * A * a1_base * a2_base))) / (2. * g);
    const double a2 =
        (t + std::sqrt(t * t - (4. * g * A * a1_base * a2_base))) / (2. * g);
    const double b1 = (a2 - ((A * (b1_base * a1_base + b2_base * a2_base)) / g)) *
                      (1. / (a2 - a1));
    const double b2 = 1. - b1;

    return this->two_term_xs(a1, a2, b1, b2, Ee, ndl);
  } else {
    return this->two_term_xs(a1_base, a2_base, b1_base, b2_base, Ee, ndl);
  }
}

std::shared_ptr<CrossSection> Material::dilution_xs(
    const std::vector<double>& dils, std::shared_ptr<NDLibrary> ndl) const {
  if (dils.size() != this->size()) {
    std::stringstream mssg;
    mssg << "The number of provided dilutions (" << dils.size()
         << ") does not match the number of nuclides (" << this->size() << ").";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  for (std::size_t i = 0; i < dils.size(); i++) {
    if (dils[i] <= 0.) {
      std::stringstream mssg;
      mssg << "The provided dilution at index " << i << " is negative.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  // Get the first nuclide, then add others to it
  std::shared_ptr<CrossSection> xsout(nullptr);

  // Add other components
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const double Ni = this->atom_density(namei);
    auto xsi = ndl->interp_xs(namei, temperature_, dils[i]);
    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  return xsout;
}

std::shared_ptr<CrossSection> Material::ring_carlvik_xs(double C, double Rpin, double Rin, double Rout, std::shared_ptr<NDLibrary> ndl) const {
  const double a1 = 0.5 * (C + 5. - std::sqrt(C * C + 34. * C + 1.));
  const double a2 = 0.5 * (C + 5. + std::sqrt(C * C + 34. * C + 1.));
  const double b1 = (a2 - (1. - C)) / (a2 - a1);
  const double b2 = 1. - b1;

  std::shared_ptr<CrossSection> xsout(nullptr);

  const double mat_pot_xs = this->potential_xs();
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const double Ni = this->atom_density(namei);
    
    auto xsi = ndl->ring_two_term_xs(namei, temperature(), a1, a2, b1, b2, mat_pot_xs, Ni, Rpin, Rin, Rout);
    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  return xsout;
}

std::shared_ptr<CrossSection> Material::two_term_xs(
    const double a1, const double a2, const double b1, const double b2,
    const double Ee, std::shared_ptr<NDLibrary> ndl) const {
  std::shared_ptr<CrossSection> xsout(nullptr);

  // Add component from each nuclide
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const double Ni = this->atom_density(namei);
    const double pot_xs = ndl->get_nuclide(namei).potential_xs;
    const double macro_pot_xs = Ni * pot_xs;
    const double bg_xs_1 = (potential_xs() - macro_pot_xs + a1 * Ee) / Ni;
    const double bg_xs_2 = (potential_xs() - macro_pot_xs + a2 * Ee) / Ni;

    auto xsi = ndl->two_term_xs(namei, temperature(), b1, b2, bg_xs_1, bg_xs_2);
    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  return xsout;
}

}  // namespace scarabee

// References
// [1] H. Koike, K. Yamaji, K. Kirimura, D. Sato, H. Matsumoto, and A. Yamamoto,
//     “Advanced resonance self-shielding method for gray resonance treatment in
//     lattice physics code GALAXY,” J. Nucl. Sci. Technol., vol. 49, no. 7,
//     pp. 725–747, 2012, doi: 10.1080/00223131.2012.693885.
//
// [2] N. Gibson, “Novel Resonance Self-Shielding Methods for Nuclear Reactor
//     Analysis,” Massachusetts Institute of Technology, 2016.
