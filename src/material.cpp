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

  components.emplace_back();
  components.back().name = name;
  components.back().fraction = frac;
}

void MaterialComposition::add_nuclide(const MaterialComponent& comp) {
  if (comp.fraction <= 0.) {
    std::stringstream mssg;
    mssg << "Component \"" << comp.name << "\" has negative or zero fraction.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  components.push_back(comp);
}

Material::Material(const MaterialComposition& comp, double temp, double density, DensityUnits du, std::shared_ptr<NDLibrary> ndl):
composition_(comp),
dilutions_(comp.components.size(), 1.E10), // Default to infinite dilution
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

  if (density < 0.) {
    auto mssg = "Material density must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (composition_.components.empty()) {
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
    for (const auto& c : composition_.components)
      frac_sum += c.fraction;
    
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
    for (auto& c : composition_.components) {
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
  for (const auto& c : composition_.components) {
    const auto& nuc = ndl->get_nuclide(c.name);
    if (nuc.fissile)
      fissile_ = true;

    if (nuc.resonant)
      resonant_ = true;

    potential_xs_ += atoms_per_bcm_ * c.fraction * nuc.potential_xs;
  }
}

void Material::set_dilution(const std::string& name, double d) {
  if (d <= 0.) {
    auto mssg = "Dilution must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < composition_.components.size(); i++) {
    if (composition_.components[i].name == name) {
      dilutions_[i] = d;
      return;
    }
  }

  std::stringstream mssg;
  mssg << "Could not find component with name \"" << name << "\".";
  spdlog::error(mssg.str());
  throw ScarabeeException(mssg.str());
}

double Material::calc_avg_molar_mass(const NDLibrary& ndl) const {
  double avg_mm = 0.;

  for (const auto& comp : composition_.components) {
    const auto& nuc = ndl.get_nuclide(comp.name);
    if (composition_.fractions == Fraction::Atoms) {
      avg_mm += comp.fraction * nuc.awr * N_MASS_AMU;
    } else {
      avg_mm += comp.fraction / (nuc.awr * N_MASS_AMU);
    }
  }

  if (composition_.fractions == Fraction::Weight)
    avg_mm = 1. / avg_mm;

  return avg_mm;
}

void Material::normalize_fractions() {
  double frac_sum = 0.;

  for (const auto& comp : composition_.components) {
    frac_sum += comp.fraction;
  }

  for (auto& comp : composition_.components) {
    comp.fraction /= frac_sum;
  }
}

}
