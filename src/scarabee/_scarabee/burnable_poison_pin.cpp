#include <assemblies/pins/burnable_poison_pin.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <memory>

namespace scarabee {

BurnablePoisonPin::BurnablePoisonPin(
    std::shared_ptr<Material> center, double center_radius,
    std::shared_ptr<Material> poison_clad, double inner_poison_clad_radius,
    std::shared_ptr<Material> gap, std::optional<double> inner_gap_radius,
    std::shared_ptr<Material> poison, double poison_radius,
    std::optional<double> outer_gap_radius, double outer_poison_clad_radius,
    double inner_moderator_radius, std::shared_ptr<Material> guide_tube_clad,
    double guide_tube_radius)
    : center_(center),
      center_radius_(center_radius),
      poison_clad_(poison_clad),
      inner_poison_clad_radius_(inner_poison_clad_radius),
      gap_(gap),
      inner_gap_radius_(0.),
      poison_(poison),
      poison_radius_(poison_radius),
      outer_gap_radius_(0.),
      outer_poison_clad_radius_(outer_poison_clad_radius),
      inner_moderator_radius_(inner_moderator_radius),
      guide_tube_clad_(guide_tube_clad),
      guide_tube_radius_(guide_tube_radius),
      condensed_xs_() {
  if (gap && inner_gap_radius && outer_gap_radius) {
    inner_gap_radius_ = inner_gap_radius.value();
    outer_gap_radius_ = outer_gap_radius.value();

    if (inner_gap_radius_ >= outer_gap_radius_) {
      auto mssg = "Inner gap radius must be > outer gap radius.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  } else if (!gap && !inner_gap_radius && !outer_gap_radius) {
    // We are good, no gap and no data for other bits
  } else {
    auto mssg = "Must provide all three gap parameters.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

Vector BurnablePoisonPin::clad_offset() const {
  return Vector(inner_moderator_radius_ +
                    0.5 * (guide_tube_radius_ - inner_moderator_radius_),
                0.0);
}

std::shared_ptr<SimplePinCell> BurnablePoisonPin::make_fuel_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  std::vector<double> radii;
  std::vector<std::shared_ptr<CrossSection>> mats;
  radii.reserve(9);
  mats.reserve(9);

  xt::xtensor<double, 1> Et{center_->potential_xs()};
  xt::xtensor<double, 1> Ea{Et(0)};
  xt::xtensor<double, 2> Es{{0.}};
  radii.push_back(center_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "BP Center"));

  Et(0) = poison_clad_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(inner_poison_clad_radius_);
  auto PC = std::make_shared<CrossSection>(Et, Ea, Es, "Poison Clad");
  mats.push_back(PC);

  std::shared_ptr<CrossSection> Gap{nullptr};
  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    radii.push_back(inner_gap_radius_);
    Gap = std::make_shared<CrossSection>(Et, Ea, Es, "Gap");
    mats.push_back(Gap);
  }

  Et(0) = poison_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(poison_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Poison"));

  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    radii.push_back(outer_gap_radius_);
    mats.push_back(Gap);
  }

  Et(0) = poison_clad_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(outer_poison_clad_radius_);
  mats.push_back(PC);

  Et(0) = moderator->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(inner_moderator_radius_);
  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Moderator");
  mats.push_back(Mod);

  Et(0) = guide_tube_clad_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(guide_tube_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Clad"));

  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<SimplePinCell> BurnablePoisonPin::make_clad_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  std::vector<double> radii;
  std::vector<std::shared_ptr<CrossSection>> mats;
  radii.reserve(9);
  mats.reserve(9);

  xt::xtensor<double, 1> Et{center_->potential_xs()};
  xt::xtensor<double, 1> Ea{Et(0)};
  xt::xtensor<double, 2> Es{{0.}};
  radii.push_back(center_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "BP Center"));

  Et(0) = poison_clad_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(inner_poison_clad_radius_);
  auto PC = std::make_shared<CrossSection>(Et, Ea, Es, "Poison Clad");
  mats.push_back(PC);

  std::shared_ptr<CrossSection> Gap{nullptr};
  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    radii.push_back(inner_gap_radius_);
    Gap = std::make_shared<CrossSection>(Et, Ea, Es, "Gap");
    mats.push_back(Gap);
  }

  Et(0) = poison_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(poison_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Poison"));

  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    radii.push_back(outer_gap_radius_);
    mats.push_back(Gap);
  }

  Et(0) = poison_clad_->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(outer_poison_clad_radius_);
  mats.push_back(PC);

  Et(0) = moderator->potential_xs();
  Ea(0) = Et(0);
  radii.push_back(inner_moderator_radius_);
  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Moderator");
  mats.push_back(Mod);

  Et(0) = 1.0E5;
  Ea(0) = Et(0);
  radii.push_back(guide_tube_radius_);
  mats.push_back(std::make_shared<CrossSection>(Et, Ea, Es, "Clad"));

  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<CylindricalCell> BurnablePoisonPin::make_cylindrical_cell(
    double pitch, std::shared_ptr<CrossSection> moderator, double buffer_radius,
    std::shared_ptr<CrossSection> buffer, std::shared_ptr<NDLibrary> ndl,
    std::optional<double> dancoff_clad, double clad_dilution,
    double poison_clad_dilution) const {
  std::vector<double> radii;
  std::vector<std::shared_ptr<CrossSection>> mats;
  radii.reserve(10);
  mats.reserve(10);

  radii.push_back(center_radius_);
  mats.push_back(
      center_->dilution_xs(std::vector<double>(center_->size(), 1.E10), ndl));
  if (mats.back()->name() == "") mats.back()->set_name("Center");

  radii.push_back(inner_poison_clad_radius_);
  auto PC = poison_clad_->dilution_xs(
      std::vector<double>(poison_clad_->size(), poison_clad_dilution), ndl);
  mats.push_back(PC);
  if (mats.back()->name() == "") mats.back()->set_name("Poison Clad");

  std::shared_ptr<CrossSection> Gap{nullptr};
  if (gap_) {
    radii.push_back(inner_gap_radius_);
    Gap = gap_->dilution_xs(std::vector<double>(gap_->size(), 1.0E10), ndl);
    mats.push_back(Gap);
    if (mats.back()->name() == "") mats.back()->set_name("Gap");
  }

  radii.push_back(poison_radius_);
  mats.push_back(
      poison_->dilution_xs(std::vector<double>(poison_->size(), 1.0E10), ndl));
  if (mats.back()->name() == "") mats.back()->set_name("Poison");

  if (gap_) {
    radii.push_back(outer_gap_radius_);
    mats.push_back(Gap);
  }

  radii.push_back(outer_poison_clad_radius_);
  mats.push_back(PC);

  radii.push_back(inner_moderator_radius_);
  mats.push_back(moderator);

  radii.push_back(guide_tube_radius_);
  if (dancoff_clad) {
    double Ee = 1. / (2. * (guide_tube_radius_ - inner_moderator_radius_));
    mats.push_back(guide_tube_clad_->roman_xs(dancoff_clad.value(), Ee, ndl));
  } else {
    mats.push_back(guide_tube_clad_->dilution_xs(
        std::vector<double>(guide_tube_clad_->size(), clad_dilution), ndl));
  }
  if (mats.back()->name() == "") mats.back()->set_name("Clad");

  radii.push_back(std::sqrt(pitch * pitch / PI));
  mats.push_back(moderator);

  // Add the buffer
  radii.push_back(buffer_radius);
  mats.push_back(buffer);

  return std::make_shared<CylindricalCell>(radii, mats);
}

std::shared_ptr<PinCell> BurnablePoisonPin::make_moc_cell(double pitch) const {
  std::vector<double> radii;
  radii.reserve(8);

  radii.push_back(center_radius_);
  radii.push_back(inner_poison_clad_radius_);
  if (gap_) radii.push_back(inner_gap_radius_);
  radii.push_back(poison_radius_);
  if (gap_) radii.push_back(outer_gap_radius_);
  radii.push_back(outer_poison_clad_radius_);
  radii.push_back(inner_moderator_radius_);
  radii.push_back(guide_tube_radius_);

  return std::make_shared<PinCell>(radii, condensed_xs_, pitch, pitch);
}

void BurnablePoisonPin::load_nuclides(std::shared_ptr<NDLibrary> ndl) const {
  center_->load_nuclides(ndl);
  poison_clad_->load_nuclides(ndl);
  poison_->load_nuclides(ndl);
  guide_tube_clad_->load_nuclides(ndl);
  if (gap_) gap_->load_nuclides(ndl);
}

}  // namespace scarabee
