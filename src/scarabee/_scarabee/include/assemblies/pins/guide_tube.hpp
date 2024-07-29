#ifndef SCARABEE_GUIDE_TUBE_H
#define SCARABEE_GUIDE_TUBE_H

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

class GuideTube {
 public:
  GuideTube(std::shared_ptr<Material> clad, double inner_radius, double outer_radius);

  Vector clad_offset() const;

  std::shared_ptr<SimplePinCell> make_fuel_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<SimplePinCell> make_clad_dancoff_cell(
      double pitch, std::shared_ptr<Material> moderator) const;

  std::shared_ptr<CylindricalCell> make_cylindrical_cell(
      double pitch, std::shared_ptr<CrossSection> moderator,
      double buffer_radius,
      std::shared_ptr<CrossSection> buffer,
      std::shared_ptr<NDLibrary> ndl,
      std::optional<double> dancoff_clad, double clad_dilution = 300.) const;

  std::shared_ptr<PinCell> make_moc_cell(double pitch) const;

  std::vector<std::shared_ptr<CrossSection>>& condensed_xs() {
    return condensed_xs_;
  }
  const std::vector<std::shared_ptr<CrossSection>>& condensed_xs() const {
    return condensed_xs_;
  }

 private:
  std::shared_ptr<Material> clad_;
  double inner_radius_;
  double outer_radius_;
  std::vector<std::shared_ptr<CrossSection>> condensed_xs_;
};

}  // namespace scarabee

#endif
