#include <assemblies/pins/fuel_pin.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>

namespace scarabee {

FuelPin::FuelPin(std::shared_ptr<Material> fuel, double fuel_radius,
                 std::shared_ptr<Material> gap,
                 std::optional<double> gap_radius,
                 std::shared_ptr<Material> clad, double clad_radius,
                 std::size_t fuel_rings, bool needs_buffer)
    : fuel_(fuel),
      fuel_radius_(fuel_radius),
      gap_(gap),
      gap_radius_(0.),
      clad_(clad),
      clad_radius_(clad_radius),
      fuel_rings_(fuel_rings),
      condensed_xs_(),
      needs_buffer_(needs_buffer) {
  if (gap_ == nullptr && gap_radius.has_value()) {
    auto mssg = "Fuel gap material is provided, but no gap radius.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (gap_ && gap_radius.has_value() == false) {
    auto mssg = "Fuel gap radius provided, but no gap material.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (gap_) {
    gap_radius_ = gap_radius.value();
  }

  if (fuel_rings_ == 0 || fuel_rings_ > 20) {
    auto mssg = "Number of fuel rings must be >= 1 and <= 20.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

Vector FuelPin::clad_offset() const {
  if (gap_) {
    return Vector(gap_radius_ + 0.5 * (clad_radius_ - gap_radius_), 0.);
  } else {
    return Vector(gap_radius_ + 0.5 * (clad_radius_ - fuel_radius_), 0.);
  }
}

std::shared_ptr<SimplePinCell> FuelPin::make_fuel_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  // We first determine all the radii
  std::vector<double> radii;
  radii.reserve(3);
  radii.push_back(fuel_radius_);
  if (gap_) radii.push_back(gap_radius_);
  radii.push_back(clad_radius_);

  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(4);

  // This returns a cell for calculating the fuel pin Dancoff factor.
  // As such, the fuel XS has infinite values.
  xt::xtensor<double, 1> Et = {1.0E5};
  xt::xtensor<double, 1> Ea = {Et(0)};
  xt::xtensor<double, 2> Es = {{0.}};
  auto Fuel = std::make_shared<CrossSection>(Et, Ea, Es, "Fuel");
  mats.push_back(Fuel);

  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    auto Gap = std::make_shared<CrossSection>(Et, Ea, Es, "Gap");
    mats.push_back(Gap);
  }

  Et(0) = clad_->potential_xs();
  Ea(0) = Et(0);
  auto Clad = std::make_shared<CrossSection>(Et, Ea, Es, "Clad");
  mats.push_back(Clad);

  Et(0) = moderator->potential_xs();
  Ea(0) = Et(0);
  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Moderator");
  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<SimplePinCell> FuelPin::make_clad_dancoff_cell(
    double pitch, std::shared_ptr<Material> moderator) const {
  // We first determine all the radii
  std::vector<double> radii;
  radii.reserve(3);
  radii.push_back(fuel_radius_);
  if (gap_) radii.push_back(gap_radius_);
  radii.push_back(clad_radius_);

  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(4);

  xt::xtensor<double, 1> Et = {fuel_->potential_xs()};
  xt::xtensor<double, 1> Ea = {Et(0)};
  xt::xtensor<double, 2> Es = {{0.}};
  auto Fuel = std::make_shared<CrossSection>(Et, Ea, Es, "Fuel");
  mats.push_back(Fuel);

  if (gap_) {
    Et(0) = gap_->potential_xs();
    Ea(0) = Et(0);
    auto Gap = std::make_shared<CrossSection>(Et, Ea, Es, "Gap");
    mats.push_back(Gap);
  }

  // This returns a cell for calculating the clad Dancoff factor.
  // As such, the clad XS has infinite values.
  Et(0) = 1.0E5;
  Ea(0) = Et(0);
  auto Clad = std::make_shared<CrossSection>(Et, Ea, Es, "Clad");
  mats.push_back(Clad);

  Et(0) = moderator->potential_xs();
  Ea(0) = Et(0);
  auto Mod = std::make_shared<CrossSection>(Et, Ea, Es, "Moderator");
  mats.push_back(Mod);

  return std::make_shared<SimplePinCell>(radii, mats, pitch, pitch);
}

std::shared_ptr<CylindricalCell> FuelPin::make_cylindrical_cell(
    double pitch, double dancoff_fuel, std::shared_ptr<CrossSection> moderator,
    std::shared_ptr<NDLibrary> ndl, std::optional<double> dancoff_clad,
    double clad_dilution) const {
  // We first determine all the radii
  std::vector<double> radii;
  radii.reserve(fuel_rings_ + 3);

  // Add the fuel radii
  if (fuel_rings_ == 1) {
    radii.push_back(fuel_radius_);
  } else {
    // We will subdivide the pellet into rings. Each ring will have the same
    // volume. Start by getting pellet volume.
    const double V = PI * fuel_radius_ * fuel_radius_;
    const double Vr = V / static_cast<double>(fuel_rings_);
    for (std::size_t r = 0; r < fuel_rings_; r++) {
      const double Rin = r > 0 ? radii.back() : 0.;
      double Rout = std::sqrt((Vr + PI * Rin * Rin) / PI);
      if (Rout > fuel_radius_) Rout = fuel_radius_;
      radii.push_back(Rout);
    }
  }

  // Add gap radius
  if (gap_) radii.push_back(gap_radius_);

  // Add clad radius
  radii.push_back(clad_radius_);

  // Add water radius
  radii.push_back(std::sqrt(pitch * pitch / PI));

  // Next, we determine all the materials.
  std::vector<std::shared_ptr<CrossSection>> mats;

  // This requires applying self shielding to the fuel and cladding
  if (fuel_rings_ == 1) {
    // For a single pellet, use the standard Carlvik two-term
    double Ee = 1.0 / (2.0 * fuel_radius_);  // Fuel escape xs
    mats.push_back(fuel_->carlvik_xs(dancoff_fuel, Ee, ndl));
    if (mats.back()->name() == "") mats.back()->set_name("Fuel");
  } else {
    // We need to apply spatial self shielding
    for (std::size_t r = 0; r < fuel_rings_; r++) {
      const double Rin = r > 0 ? radii[r - 1] : 0.;
      const double Rout = radii[r];
      mats.push_back(
          fuel_->ring_carlvik_xs(dancoff_fuel, fuel_radius_, Rin, Rout, ndl));
      if (mats.back()->name() == "") mats.back()->set_name("Fuel");
    }
  }

  // Next, add the gap (if present)
  if (gap_) {
    std::vector<double> dilutions(gap_->size(), 1.0E10);
    mats.push_back(gap_->dilution_xs(dilutions, ndl));
    if (mats.back()->name() == "") mats.back()->set_name("Gap");
  }

  // Add the cladding
  if (dancoff_clad) {
    double Ee = 0.;
    if (gap_) {
      Ee = 1. / (2. * (clad_radius_ - gap_radius_));
    } else {
      Ee = 1. / (2. * (clad_radius_ - fuel_radius_));
    }
    mats.push_back(clad_->roman_xs(dancoff_clad.value(), Ee, ndl));
  } else {
    std::vector<double> dilutions(clad_->size(), clad_dilution);
    mats.push_back(clad_->dilution_xs(dilutions, ndl));
  }
  if (mats.back()->name() == "") mats.back()->set_name("Clad");

  // Finally, add moderator
  mats.push_back(moderator);

  return std::make_shared<CylindricalCell>(radii, mats);
}

std::shared_ptr<CylindricalCell> FuelPin::make_cylindrical_cell(
      double pitch, double buffer_radius, std::shared_ptr<CrossSection> buffer,
      double dancoff_fuel, std::shared_ptr<CrossSection> moderator,
      std::shared_ptr<NDLibrary> ndl, std::optional<double> dancoff_clad,
      double clad_dilution) const {
  // We first determine all the radii
  std::vector<double> radii;
  radii.reserve(fuel_rings_ + 4);

  // Add the fuel radii
  if (fuel_rings_ == 1) {
    radii.push_back(fuel_radius_);
  } else {
    // We will subdivide the pellet into rings. Each ring will have the same
    // volume. Start by getting pellet volume.
    const double V = PI * fuel_radius_ * fuel_radius_;
    const double Vr = V / static_cast<double>(fuel_rings_);
    for (std::size_t r = 0; r < fuel_rings_; r++) {
      const double Rin = r > 0 ? radii.back() : 0.;
      double Rout = std::sqrt((Vr + PI * Rin * Rin) / PI);
      if (Rout > fuel_radius_) Rout = fuel_radius_;
      radii.push_back(Rout);
    }
  }

  // Add gap radius
  if (gap_) radii.push_back(gap_radius_);

  // Add clad radius
  radii.push_back(clad_radius_);

  // Add water radius
  radii.push_back(std::sqrt(pitch * pitch / PI));

  // Add buffer radius
  radii.push_back(buffer_radius);

  // Next, we determine all the materials.
  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(fuel_rings_ + 4);

  // This requires applying self shielding to the fuel and cladding
  if (fuel_rings_ == 1) {
    // For a single pellet, use the standard Carlvik two-term
    double Ee = 1.0 / (2.0 * fuel_radius_);  // Fuel escape xs
    mats.push_back(fuel_->carlvik_xs(dancoff_fuel, Ee, ndl));
    if (mats.back()->name() == "") mats.back()->set_name("Fuel");
  } else {
    // We need to apply spatial self shielding
    for (std::size_t r = 0; r < fuel_rings_; r++) {
      const double Rin = r > 0 ? radii[r - 1] : 0.;
      const double Rout = radii[r];
      mats.push_back(
          fuel_->ring_carlvik_xs(dancoff_fuel, fuel_radius_, Rin, Rout, ndl));
      if (mats.back()->name() == "") mats.back()->set_name("Fuel");
    }
  }

  // Next, add the gap (if present)
  if (gap_) {
    std::vector<double> dilutions(gap_->size(), 1.0E10);
    mats.push_back(gap_->dilution_xs(dilutions, ndl));
    if (mats.back()->name() == "") mats.back()->set_name("Gap");
  }

  // Add the cladding
  if (dancoff_clad) {
    double Ee = 0.;
    if (gap_) {
      Ee = 1. / (2. * (clad_radius_ - gap_radius_));
    } else {
      Ee = 1. / (2. * (clad_radius_ - fuel_radius_));
    }
    mats.push_back(clad_->roman_xs(dancoff_clad.value(), Ee, ndl));
  } else {
    std::vector<double> dilutions(clad_->size(), clad_dilution);
    mats.push_back(clad_->dilution_xs(dilutions, ndl));
  }
  if (mats.back()->name() == "") mats.back()->set_name("Clad");

  // Add moderator
  mats.push_back(moderator);

  // Finally, add the buffer
  mats.push_back(buffer);

  return std::make_shared<CylindricalCell>(radii, mats);     
}

std::shared_ptr<PinCell> FuelPin::make_moc_cell(double pitch) const {
  // We first determine all the radii
  std::vector<double> radii;
  radii.reserve(fuel_rings_ + 4);

  // Add the fuel radii
  if (fuel_rings_ == 1) {
    radii.push_back(fuel_radius_);
  } else {
    // We will subdivide the pellet into rings. Each ring will have the same
    // volume. Start by getting pellet volume.
    const double V = PI * fuel_radius_ * fuel_radius_;
    const double Vr = V / static_cast<double>(fuel_rings_);
    for (std::size_t r = 0; r < fuel_rings_; r++) {
      const double Rin = r > 0 ? radii.back() : 0.;
      double Rout = std::sqrt((Vr + PI * Rin * Rin) / PI);
      if (Rout > fuel_radius_) Rout = fuel_radius_;
      radii.push_back(Rout);
    }
  }

  // Add gap radius
  if (gap_) radii.push_back(gap_radius_);

  // Add clad radius
  radii.push_back(clad_radius_);

  // Add intermediate water radius
  const double mod_width = 0.5 * pitch - radii.back();
  radii.push_back(radii.back() + 0.8 * mod_width);

  std::vector<std::shared_ptr<CrossSection>> mats;
  mats.reserve(condensed_xs_.size() + 1);
  mats = condensed_xs_;
  mats.push_back(condensed_xs_.back());

  return std::make_shared<PinCell>(radii, mats, pitch, pitch);
}

void FuelPin::load_nuclides(std::shared_ptr<NDLibrary> ndl) const {
  fuel_->load_nuclides(ndl);
  clad_->load_nuclides(ndl);
  if (gap_) gap_->load_nuclides(ndl);
}

}  // namespace scarabee
