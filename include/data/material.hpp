#ifndef SCARABEE_MATERIAL_H
#define SCARABEE_MATERIAL_H

#include <string>
#include <memory>
#include <vector>

namespace scarabee {

struct MaterialComponent {
  std::string name;
  double fraction;
};

enum class Fraction {Atoms, Weight};

struct MaterialComposition {
  std::vector<MaterialComponent> components;
  Fraction fractions;

  void add_nuclide(const std::string& name, double frac);
  void add_nuclide(const MaterialComponent& comp);
};

enum class DensityUnits {g_cm3, a_bcm, sum};

// Pre-declared for converting fractions
class NDLibrary;

class Material {
  public:
    Material(const MaterialComposition& comp, double temp, std::shared_ptr<NDLibrary> ndl);
    Material(const MaterialComposition& comp, double temp, double density, DensityUnits du, std::shared_ptr<NDLibrary> ndl);

    const MaterialComposition& composition() const { return composition_; }

    bool has_component(const std::string& name) const;
    double atom_density(const std::string& name) const;

    double atoms_per_bcm() const { return atoms_per_bcm_; }
    double grams_per_cm3() const { return grams_per_cm3_; }
    double potential_xs() const { return potential_xs_; }
    double temperature() const { return temperature_; }
    double average_molar_mass() const { return average_molar_mass_; }

    bool fissile() const { return fissile_; }
    bool resonant() const { return resonant_; }

  private:
    MaterialComposition composition_;
    double temperature_;
    double average_molar_mass_;
    double atoms_per_bcm_;
    double grams_per_cm3_;
    double potential_xs_;
    bool fissile_;
    bool resonant_;

    double calc_avg_molar_mass(const NDLibrary& ndl) const;
    void normalize_fractions();
};

}

#endif