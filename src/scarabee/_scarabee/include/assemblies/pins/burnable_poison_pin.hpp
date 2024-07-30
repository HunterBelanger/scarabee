#ifndef SCARABEE_BURNABLE_POISON_PIN_H
#define SCARABEE_BURNABLE_POISON_PIN_H

#include <data/material.hpp>
#include <data/nd_library.hpp>
#include <cross_section.hpp>
#include <moc/vector.hpp>
#include <moc/cell.hpp>
#include <moc/simple_pin_cell.hpp>
#include <moc/pin_cell.hpp>
#include <cylindrical_cell.hpp>

#include <memory>
#include <optional>
#include <vector>

namespace scarabee {

class BurnablePoisonPin {
 public:
  BurnablePoisonPin(
      std::shared_ptr<Material> center, double center_radius,
      std::shared_ptr<Material> poison_clad, double inner_poison_clad_radius,
      std::shared_ptr<Material> gap, std::optional<double> inner_gap_radius,
      std::shared_ptr<Material> poison, double poison_radius,
      std::optional<double> outer_gap_radius, double outer_poison_clad_radius,
      double inner_moderator_radius, std::shared_ptr<Material> guide_tube_clad,
      double guide_tube_radius);

  Vector clad_offset() const;

  std::shared_ptr<SimplePinCell> make_fuel_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<SimplePinCell> make_clad_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<CylindricalCell> make_cylindrical_cell(
      double pitch, std::shared_ptr<CrossSection> moderator,
      double buffer_radius, std::shared_ptr<CrossSection> buffer,
      std::shared_ptr<NDLibrary> ndl, std::optional<double> dancoff_clad,
      double clad_dilution = 300., double poison_clad_dilution = 300) const;

  std::shared_ptr<PinCell> make_moc_cell(double pitch) const;

  std::vector<std::shared_ptr<CrossSection>>& condensed_xs() {
    return condensed_xs_;
  }

  const std::vector<std::shared_ptr<CrossSection>>& condensed_xs() const {
    return condensed_xs_;
  }

  void load_nuclides(std::shared_ptr<NDLibrary> ndl) const;

 private:
  std::shared_ptr<Material> center_;
  double center_radius_;
  std::shared_ptr<Material> poison_clad_;
  double inner_poison_clad_radius_;
  std::shared_ptr<Material> gap_;
  double inner_gap_radius_;
  std::shared_ptr<Material> poison_;
  double poison_radius_;
  double outer_gap_radius_;
  double outer_poison_clad_radius_;
  double inner_moderator_radius_;
  std::shared_ptr<Material> guide_tube_clad_;
  double guide_tube_radius_;
  std::vector<std::shared_ptr<CrossSection>> condensed_xs_;
};

}  // namespace scarabee

#endif
