#ifndef SCARABEE_MATERIAL_H
#define SCARABEE_MATERIAL_H

#include <data/cross_section.hpp>

#include <string>
#include <memory>
#include <optional>
#include <vector>

namespace scarabee {

struct Nuclide {
  std::string name;
  double fraction;
};

enum class Fraction { Atoms, Weight };

struct MaterialComposition {
  std::vector<Nuclide> nuclides;
  Fraction fractions;
  std::string name;

  MaterialComposition(Fraction f = Fraction::Atoms, const std::string& name = "");

  void add_element(const std::string& name, double frac);
  void add_leu(double enrichment, double frac);
  void add_nuclide(const std::string& name, double frac);
  void add_nuclide(const Nuclide& nuc);
};

enum class DensityUnits { g_cm3, a_bcm, sum };

// Pre-declared for converting fractions
class NDLibrary;

class Material {
 public:
  Material(const MaterialComposition& comp, double temp,
           std::shared_ptr<NDLibrary> ndl);
  Material(const MaterialComposition& comp, double temp, double density,
           DensityUnits du, std::shared_ptr<NDLibrary> ndl);

  const MaterialComposition& composition() const { return composition_; }

  std::size_t size() const { return composition_.nuclides.size(); }

  const std::string& name() const { return name_; }
  void set_name(const std::string& new_name) { name_ = new_name; }

  std::size_t max_legendre_order() const { return max_l_; }
  void set_max_legendre_order(std::size_t max_l) { max_l_ = max_l; }

  bool has_component(const std::string& name) const;
  double atom_density(const std::string& name) const;

  double atoms_per_bcm() const { return atoms_per_bcm_; }
  double grams_per_cm3() const { return grams_per_cm3_; }
  double potential_xs() const { return potential_xs_; }
  double lambda_potential_xs() const { return lambda_potential_xs_; }
  double temperature() const { return temperature_; }
  double average_molar_mass() const { return average_molar_mass_; }

  bool fissile() const { return fissile_; }
  bool resonant() const { return resonant_; }

  std::shared_ptr<CrossSection> carlvik_xs(
      double C, double Ee, std::shared_ptr<NDLibrary> ndl,
      std::optional<std::size_t> max_l = std::nullopt) const;

  std::shared_ptr<CrossSection> roman_xs(
      double C, double Ee, std::shared_ptr<NDLibrary> ndl,
      std::optional<std::size_t> max_l = std::nullopt) const;

  std::shared_ptr<CrossSection> dilution_xs(
      const std::vector<double>& dils, std::shared_ptr<NDLibrary> ndl,
      std::optional<std::size_t> max_l = std::nullopt) const;

  std::shared_ptr<CrossSection> ring_carlvik_xs(
      double C, double Rfuel, double Rin, double Rout,
      std::shared_ptr<NDLibrary> ndl,
      std::optional<std::size_t> max_l = std::nullopt) const;

  void load_nuclides(std::shared_ptr<NDLibrary> ndl) const;

 private:
  MaterialComposition composition_;
  std::string name_;
  double temperature_;
  double average_molar_mass_;
  double atoms_per_bcm_;
  double grams_per_cm3_;
  double potential_xs_;
  double lambda_potential_xs_;
  std::size_t max_l_{1};
  bool fissile_;
  bool resonant_;

  double calc_avg_molar_mass(const NDLibrary& ndl) const;
  void normalize_fractions();

  std::shared_ptr<CrossSection> two_term_xs(const double a1, const double a2,
                                            const double b1, const double b2,
                                            const double Ee,
                                            std::shared_ptr<NDLibrary> ndl,
                                            std::size_t max_l) const;
};

}  // namespace scarabee

#endif
