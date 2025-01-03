#include <data/material.hpp>
#include <data/nd_library.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <algorithm>
#include <map>
#include <sstream>

namespace scarabee {

std::string nuclide_name_to_simple_name(const std::string& name) {
  // Must remove any _ from nuclide name (for TSLs like H1_H2O)
  std::string simp_name = name;
  auto loc_undr_scr = simp_name.find('_');
  if (loc_undr_scr != std::string::npos) {
    simp_name.resize(loc_undr_scr);
  }
  return simp_name;
}

MaterialComposition::MaterialComposition(Fraction f, const std::string& name)
    : nuclides(), fractions(f), name(name) {}

void MaterialComposition::add_element(const std::string& name, double frac) {
  if (frac <= 0.) {
    std::stringstream mssg;
    mssg << "Element \"" << name << "\" given negative or zero fraction.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // First, try to get list of natural isotopes
  if (ELEMENT_ISOTOPES.find(name) == ELEMENT_ISOTOPES.end()) {
    std::stringstream mssg;
    mssg << "Could not find element by name of \"" << name << "\".";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }
  const auto& isotopes = ELEMENT_ISOTOPES.at(name);

  if (fractions == Fraction::Atoms) {
    // Just add all isotopes if we are using atom fractions
    for (const auto& iso : isotopes) {
      this->add_nuclide(iso, NATURAL_ABUNDANCES.at(iso) * frac);
    }
  } else {
    // Must convert to weight fractions

    // First, get atomic mass of element
    double elem_atomic_mass = 0.;
    for (const auto& iso : isotopes) {
      elem_atomic_mass += ISOTOPE_MASSES.at(iso) * NATURAL_ABUNDANCES.at(iso);
    }

    // We can now get a mass abundance for each isotopes
    std::vector<double> mass_abundances(isotopes.size(), 0.);
    for (std::size_t i = 0; i < isotopes.size(); i++) {
      const auto& iso = isotopes[i];
      mass_abundances[i] = NATURAL_ABUNDANCES.at(iso) * ISOTOPE_MASSES.at(iso) /
                           elem_atomic_mass;
    }

    // Normalize the abundances
    const double mass_abundance_sum =
        std::accumulate(mass_abundances.begin(), mass_abundances.end(), 0.);

    for (auto& ma : mass_abundances) ma /= mass_abundance_sum;

    // Add all isotopes
    for (std::size_t i = 0; i < isotopes.size(); i++) {
      const auto& iso = isotopes[i];
      this->add_nuclide(iso, frac * mass_abundances[i]);
    }
  }
}

void MaterialComposition::add_leu(double enrichment, double frac) {
  // This uses the same methods as in OpenMC, and is only valid for U235
  // enrichments <= 5% by weight. Their correlations come from
  // ORNL/CSD/TM-244 found at https://doi.org/10.2172/5561567.
  if (enrichment < 0.) {
    const auto mssg = "LEU enrichment is less than zero.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (enrichment > 5.) {
    const auto mssg =
        "LEU enrichment is greater than 5%. Method is only valid for "
        "enrichments <= 5%.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (frac <= 0.) {
    const auto mssg = "LEU given negative or zero fraction.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const std::array<std::string, 4> isotopes{"U234", "U235", "U236", "U238"};

  // These are initially mass abundances !
  std::array<double, 4> abundances{0.0089 * enrichment, enrichment,
                                   0.0046 * enrichment,
                                   100.0 - 1.0135 * enrichment};
  const double mass_abundance_sum =
      std::accumulate(abundances.begin(), abundances.end(), 0.);
  for (auto& a : abundances) a /= mass_abundance_sum;

  if (this->fractions == Fraction::Weight) {
    for (std::size_t i = 0; i < 4; i++) {
      this->add_nuclide(isotopes[i], abundances[i] * frac);
    }
  } else {
    // Must convert mass to number abundances
    for (std::size_t i = 0; i < 4; i++) {
      abundances[i] /= ISOTOPE_MASSES.at(isotopes[i]);
    }

    const double sum_abundances =
        std::accumulate(abundances.begin(), abundances.end(), 0.);
    for (auto& a : abundances) a /= sum_abundances;

    // Add all isotopes
    for (std::size_t i = 0; i < 4; i++) {
      this->add_nuclide(isotopes[i], abundances[i] * frac);
    }
  }
}

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
      name_(comp.name),
      temperature_(temp),
      atoms_per_bcm_(-1.),
      grams_per_cm3_(-1.),
      potential_xs_(0.),
      lambda_potential_xs_(0.),
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
    lambda_potential_xs_ +=
        atoms_per_bcm_ * c.fraction * nuc.ir_lambda * nuc.potential_xs;

    if (nuc.fissile) fissile_ = true;

    if (nuc.resonant) resonant_ = true;
  }
}

Material::Material(const MaterialComposition& comp, double temp, double density,
                   DensityUnits du, std::shared_ptr<NDLibrary> ndl)
    : composition_(comp),
      name_(comp.name),
      temperature_(temp),
      atoms_per_bcm_(-1.),
      grams_per_cm3_(-1.),
      potential_xs_(0.),
      lambda_potential_xs_(0.),
      fissile_(false),
      resonant_(false) {
  // Make sure quantities are positive/valid
  if (temp <= 0.) {
    auto mssg = "Material temperature must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (density <= 0. && du != DensityUnits::sum) {
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
    lambda_potential_xs_ +=
        atoms_per_bcm_ * c.fraction * nuc.ir_lambda * nuc.potential_xs;

    if (nuc.fissile) fissile_ = true;

    if (nuc.resonant) resonant_ = true;
  }
}

void Material::set_temperature(double T) {
  if (T <= 0.) {
    auto mssg = "Material temperature must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  temperature_ = T;
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
    double C, double Ee, std::shared_ptr<NDLibrary> ndl,
    std::optional<std::size_t> max_l) const {
  if (max_l.has_value() == false) max_l = max_l_;

  // This implementation is based on the methods outlined by Koike and Gibson
  // in his PhD thesis [1,2]. We start by computing the coefficients for the
  // two-term rational approximation, modified according to the Dancoff
  // correction factor, C.
  const double a1 = 0.5 * (C + 5. - std::sqrt(C * C + 34. * C + 1.));
  const double a2 = 0.5 * (C + 5. + std::sqrt(C * C + 34. * C + 1.));
  const double b1 = (a2 - (1. - C)) / (a2 - a1);
  const double b2 = 1. - b1;

  return this->two_term_xs(a1, a2, b1, b2, Ee, ndl, *max_l);
}

std::shared_ptr<CrossSection> Material::roman_xs(
    double C, double Ee, std::shared_ptr<NDLibrary> ndl,
    std::optional<std::size_t> max_l) const {
  if (max_l.has_value() == false) max_l = max_l_;
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
    const double b1 =
        (a2 - ((A * (b1_base * a1_base + b2_base * a2_base)) / g)) *
        (1. / (a2 - a1));
    const double b2 = 1. - b1;

    return this->two_term_xs(a1, a2, b1, b2, Ee, ndl, *max_l);
  } else {
    return this->two_term_xs(a1_base, a2_base, b1_base, b2_base, Ee, ndl,
                             *max_l);
  }
}

std::shared_ptr<CrossSection> Material::dilution_xs(
    const std::vector<double>& dils, std::shared_ptr<NDLibrary> ndl,
    std::optional<std::size_t> max_l) const {
  if (max_l.has_value() == false) max_l = max_l_;

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
    auto xsi = ndl->interp_xs(namei, temperature_, dils[i], *max_l);
    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  // Give the xs the name of the material
  xsout->set_name(this->name_);

  return xsout;
}

std::shared_ptr<CrossSection> Material::ring_carlvik_xs(
    double C, double Rfuel, double Rin, double Rout,
    std::shared_ptr<NDLibrary> ndl, std::optional<std::size_t> max_l) const {
  if (max_l.has_value() == false) max_l = max_l_;

  if (Rin >= Rout) {
    auto mssg = "Rin must be < Rout.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Rout > Rfuel) {
    auto mssg = "Rout must be < Rfuel.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const double a1 = 0.5 * (C + 5. - std::sqrt(C * C + 34. * C + 1.));
  const double a2 = 0.5 * (C + 5. + std::sqrt(C * C + 34. * C + 1.));
  const double b1 = (a2 - (1. - C)) / (a2 - a1);
  const double b2 = 1. - b1;

  std::shared_ptr<CrossSection> xsout(nullptr);

  const double mat_pot_xs = this->lambda_potential_xs();
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const auto& nuclide = ndl->get_nuclide(namei);
    const double Ni = this->atom_density(namei);

    std::shared_ptr<CrossSection> xsi{nullptr};

    if (nuclide.resonant) {
      xsi = ndl->ring_two_term_xs(namei, temperature(), a1, a2, b1, b2,
                                  mat_pot_xs, Ni, Rfuel, Rin, Rout, *max_l);
    } else {
      const double dil = (mat_pot_xs - Ni * nuclide.potential_xs) / Ni;
      xsi = ndl->interp_xs(namei, temperature_, dil, *max_l);
    }

    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  // Give the xs the name of the material
  xsout->set_name(this->name_);

  return xsout;
}

std::shared_ptr<CrossSection> Material::two_term_xs(
    const double a1, const double a2, const double b1, const double b2,
    const double Ee, std::shared_ptr<NDLibrary> ndl, std::size_t max_l) const {
  std::shared_ptr<CrossSection> xsout(nullptr);

  const double mat_pot_xs = this->lambda_potential_xs();

  // Add component from each nuclide
  for (std::size_t i = 0; i < composition_.nuclides.size(); i++) {
    const std::string& namei = composition_.nuclides[i].name;
    const auto& nuclide = ndl->get_nuclide(namei);
    const double Ni = this->atom_density(namei);
    const double pot_xs = nuclide.ir_lambda * nuclide.potential_xs;
    const double macro_pot_xs = Ni * pot_xs;
    const double bg_xs_1 = (mat_pot_xs - macro_pot_xs + a1 * Ee) / Ni;
    const double bg_xs_2 = (mat_pot_xs - macro_pot_xs + a2 * Ee) / Ni;

    std::shared_ptr<CrossSection> xsi{nullptr};

    if (nuclide.resonant) {
      xsi = ndl->two_term_xs(namei, temperature(), b1, b2, bg_xs_1, bg_xs_2,
                             max_l);
    } else {
      const double dil = (mat_pot_xs - Ni * nuclide.potential_xs) / Ni;
      xsi = ndl->interp_xs(namei, temperature_, dil, max_l);
    }

    *xsi *= Ni;

    if (xsout) {
      *xsout += *xsi;
    } else {
      xsout = xsi;
    }
  }

  // Give the xs the name of the material
  xsout->set_name(this->name_);

  return xsout;
}

void Material::load_nuclides(std::shared_ptr<NDLibrary> ndl) const {
  for (const auto& nuc : composition_.nuclides) {
    auto& nuc_handle = ndl->get_nuclide(nuc.name);
    nuc_handle.load_xs_from_hdf5(*ndl, max_l_);
  }
}

std::shared_ptr<Material> mix_materials(
    const std::vector<std::shared_ptr<Material>>& mats,
    std::vector<double> fracs, Fraction f, std::shared_ptr<NDLibrary> ndl) {
  /* This method was directly taken from OpenMC's Python API. */

  // Make sure all fractions are positive
  for (const auto& v : fracs) {
    if (v <= 0.) {
      auto mssg = "All fractions must be > 0.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Make sure all materials are defined
  for (const auto& m : mats) {
    if (m == nullptr) {
      auto mssg = "No material can be None.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  if (fracs.size() != mats.size()) {
    auto mssg = "Materials and fractions lists must have same length.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (fracs.size() == 0) {
    auto mssg = "Materials and fractions lists must have a length > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ndl == nullptr) {
    auto mssg = "Provided NDLibrary cannot be None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // First, we normalize the fractions
  double norm = std::accumulate(fracs.begin(), fracs.end(), 0.);
  for (auto& v : fracs) v /= norm;

  std::vector<double> wgts(mats.size(), 0.);

  for (std::size_t m = 0; m < mats.size(); m++) {
    if (f == Fraction::Atoms) {
      wgts[m] =
          fracs[m] * mats[m]->average_molar_mass() / mats[m]->grams_per_cm3();
    } else {
      // Weight
      wgts[m] = fracs[m] / mats[m]->grams_per_cm3();
    }
  }

  // Normalize weights
  norm = std::accumulate(wgts.begin(), wgts.end(), 0.);
  for (auto& v : wgts) v /= norm;

  // Compute average temperature
  double avg_T = 0.;
  for (std::size_t m = 0; m < mats.size(); m++) {
    avg_T += wgts[m] * mats[m]->temperature();
  }

  std::map<std::string, double> atoms_per_cc;
  std::map<std::string, double> mass_per_cc;
  for (std::size_t m = 0; m < mats.size(); m++) {
    const auto& mat = mats[m];
    const double wgt = wgts[m];

    for (std::size_t n = 0; n < mat->composition().nuclides.size(); n++) {
      const auto& nuc_info = mat->composition().nuclides[n];
      const std::string nuc_name = nuc_info.name;
      const std::string simp_name = nuclide_name_to_simple_name(nuc_name);
      const auto& nuc_frac = nuc_info.fraction;
      const double atms_per_bcm = mat->atom_density(nuc_name);
      const double atms_per_cc = wgt * 1.E24 * atms_per_bcm;
      atoms_per_cc[nuc_name] += atms_per_cc;
      mass_per_cc[nuc_name] +=
          atms_per_cc * ISOTOPE_MASSES.at(simp_name) / (N_AVAGADRO * 1.E24);
    }
  }

  double tot_atoms_per_cc = 0.;
  for (const auto& nuc : atoms_per_cc) tot_atoms_per_cc += nuc.second;
  double density = 0.;
  for (const auto& nuc : mass_per_cc) density += nuc.second;

  // Create material name
  std::string new_name("");
  for (std::size_t m = 0; m < mats.size(); m++) {
    new_name.append(mats[m]->name());

    if (m != mats.size() - 1) {
      new_name.append("-");
    }
  }

  // Make new material
  MaterialComposition new_mat_comp(Fraction::Atoms, new_name);
  for (const auto& name_N_pair : atoms_per_cc) {
    new_mat_comp.add_nuclide(name_N_pair.first,
                             name_N_pair.second / tot_atoms_per_cc);
  }

  return std::make_shared<Material>(new_mat_comp, avg_T, density,
                                    DensityUnits::g_cm3, ndl);
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
