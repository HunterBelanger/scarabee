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

std::shared_ptr<CrossSection> Material::build_xs(
    double C, double Ee, std::shared_ptr<NDLibrary> ndl) const {
  // Get the first nuclide, then add others to it
  const std::string& name0 = composition_.nuclides[0].name;
  const double N0 = this->atom_density(name0);
  std::shared_ptr<CrossSection> xsout =
      ndl->carlvik_two_term(name0, potential_xs_, temperature_, N0, C, Ee);

  // Add other components
  for (std::size_t i = 1; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const double Ni = this->atom_density(namei);
    auto xsi =
        ndl->carlvik_two_term(namei, potential_xs_, temperature_, Ni, C, Ee);
    *xsout += *xsi;
  }

  return xsout;
}

std::shared_ptr<CrossSection> Material::build_xs(
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
  const std::string& name0 = composition_.nuclides[0].name;
  const double N0 = this->atom_density(name0);
  std::shared_ptr<CrossSection> xsout =
      ndl->interp_nuclide_xs(name0, temperature_, dils[0]);
  *xsout *= N0;

  // Add other components
  for (std::size_t i = 1; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const double Ni = this->atom_density(namei);
    auto xsi = ndl->interp_nuclide_xs(namei, temperature_, dils[i]);
    *xsi *= Ni;
    *xsout += *xsi;
  }

  return xsout;
}

}  // namespace scarabee
