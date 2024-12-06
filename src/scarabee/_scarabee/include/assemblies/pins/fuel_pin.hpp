#ifndef SCARABEE_FUEL_PIN_H
#define SCARABEE_FUEL_PIN_H

#include <data/material.hpp>
#include <data/nd_library.hpp>
#include <data/cross_section.hpp>
#include <moc/vector.hpp>
#include <moc/cell.hpp>
#include <moc/simple_pin_cell.hpp>
#include <moc/pin_cell.hpp>
#include <cylindrical_cell.hpp>

#include <memory>
#include <optional>
#include <vector>

namespace scarabee {

class FuelPin {
 public:
  FuelPin(std::shared_ptr<Material> fuel, double fuel_radius,
          std::shared_ptr<Material> gap, std::optional<double> gap_radius,
          std::shared_ptr<Material> clad, double clad_radius,
          std::size_t fuel_rings = 1);

  Vector clad_offset() const;

  std::shared_ptr<SimplePinCell> make_fuel_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<SimplePinCell> make_clad_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<CylindricalCell> make_cylindrical_cell(
      double pitch, double dancoff_fuel,
      std::shared_ptr<CrossSection> moderator, std::shared_ptr<NDLibrary> ndl,
      std::optional<double> dancoff_clad, double clad_dilution = 300.) const;

  std::shared_ptr<PinCell> make_moc_cell(double pitch) const;

  std::vector<std::shared_ptr<CrossSection>>& condensed_xs() {
    return condensed_xs_;
  }
  const std::vector<std::shared_ptr<CrossSection>>& condensed_xs() const {
    return condensed_xs_;
  }

  void load_nuclides(std::shared_ptr<NDLibrary> ndl) const;

 private:
  std::shared_ptr<Material> fuel_;
  double fuel_radius_;
  std::shared_ptr<Material> gap_;
  double gap_radius_;
  std::shared_ptr<Material> clad_;
  double clad_radius_;
  std::size_t fuel_rings_;
  std::vector<std::shared_ptr<CrossSection>> condensed_xs_;
};

}  // namespace scarabee

#endif
