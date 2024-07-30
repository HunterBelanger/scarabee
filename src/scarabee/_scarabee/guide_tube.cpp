#include <assemblies/pins/guide_tube.hpp>
#include <memory>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include "cross_section.hpp"
#include "cylindrical_flux_solver.hpp"
#include "spdlog/spdlog.h"

namespace scarabee {

GuideTube::GuideTube(std::shared_ptr<Material> clad, double inner_radius,
                     double outer_radius)
    : clad_(clad),
      inner_radius_(inner_radius),
      outer_radius_(outer_radius),
      condensed_xs_() {
  if (outer_radius_ <= inner_radius_) {
    auto mssg = "Outer radius must be > inner radius.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

Vector GuideTube::clad_offset() const {
  return Vector(0.5 * (inner_radius_ + outer_radius_), 0.);
}

std::shared_ptr<SimplePinCell> GuideTube::make_fuel_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  std::vector<double> radii{inner_radius_, outer_radius_};

  std::vector<std::shared_ptr<CrossSection>> mats;

  xt::xtensor<double, 1> Et{moderator->potential_xs()};
  xt::xtensor<double, 1> Ea{Et(0)};
  xt::xtensor<double, 2> Es{{0.}};

  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Mod");
  mats.push_back(Mod);

  Et(0) = clad_->potential_xs();
  Ea(0) = Et(0);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Clad"));

  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<SimplePinCell> GuideTube::make_clad_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  std::vector<double> radii{inner_radius_, outer_radius_};

  std::vector<std::shared_ptr<CrossSection>> mats;

  xt::xtensor<double, 1> Et{moderator->potential_xs()};
  xt::xtensor<double, 1> Ea{Et(0)};
  xt::xtensor<double, 2> Es{{0.}};

  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Mod");
  mats.push_back(Mod);

  Et(0) = 1.0E5;
  Ea(0) = Et(0);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Clad"));

  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<CylindricalCell> GuideTube::make_cylindrical_cell(
    double pitch, std::shared_ptr<CrossSection> moderator, double buffer_radius,
    std::shared_ptr<CrossSection> buffer, std::shared_ptr<NDLibrary> ndl,
    std::optional<double> dancoff_clad, double clad_dilution) const {
  std::vector<double> radii{inner_radius_, outer_radius_,
                            std::sqrt(pitch * pitch / PI), buffer_radius};

  if (radii[2] >= radii[3]) {
    auto mssg = "Buffer radius is smaller than the radius of the cell.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Next, create materials list, applying self shielding to cladding
  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(4);

  mats.push_back(moderator);

  // Add cladding
  if (dancoff_clad) {
    const double Ee = 1. / (2. * (outer_radius_ - inner_radius_));
    mats.push_back(clad_->roman_xs(dancoff_clad.value(), Ee, ndl));
  } else {
    std::vector<double> dilutions(clad_->size(), clad_dilution);
    mats.push_back(clad_->dilution_xs(dilutions, ndl));
  }

  // Add outer moderator
  mats.push_back(moderator);

  // Add the buffer
  mats.push_back(buffer);

  return std::make_shared<CylindricalCell>(radii, mats);
}

std::shared_ptr<PinCell> GuideTube::make_moc_cell(double pitch) const {
  const double r_inner_mod = std::sqrt(0.5 * inner_radius_ * inner_radius_);
  std::vector<double> radii{r_inner_mod, inner_radius_, outer_radius_};

  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(condensed_xs_.size() + 1);
  mats.push_back(condensed_xs_[0]);
  for (const auto& xs : condensed_xs_) mats.push_back(xs);

  return std::make_shared<PinCell>(radii, mats, pitch, pitch);
}

void GuideTube::load_nuclides(std::shared_ptr<NDLibrary> ndl) const {
  clad_->load_nuclides(ndl);
}

}  // namespace scarabee
