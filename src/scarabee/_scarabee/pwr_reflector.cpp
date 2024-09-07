#include <assemblies/pwr_reflector.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/criticality_spectrum.hpp>
#include <moc/moc_plotter.hpp>
#include <moc/empty_cell.hpp>
#include <diffusion/diffusion_data.hpp>

#include <xtensor/xio.hpp>

namespace scarabee {

PWRReflector::PWRReflector(double pitch, std::shared_ptr<Material> moderator,
                           std::pair<std::size_t, std::size_t> shape,
                           double gap_width, double baffle_width,
                           std::shared_ptr<Material> baffle,
                           std::shared_ptr<NDLibrary> ndl)
    : pitch_(pitch),
      shape_(shape),
      ndl_(ndl),
      gap_width_(gap_width),
      baffle_width_(baffle_width),
      baffle_(baffle) {
  this->set_moderator(moderator);

  if (pitch_ <= 0.) {
    auto mssg = "Pitch must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (gap_width_ < 0.) {
    auto mssg = "Reflector gap width must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (baffle_width_ < 0.) {
    auto mssg = "Reflector baffle width must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (baffle_ == nullptr && baffle_width_ != 0.) {
    auto mssg = "Reflector baffle material is None but baffle width is > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (baffle_width_ == 0. && baffle_ != nullptr) {
    auto mssg = "Reflector baffle width is 0, but baffle material is not None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const double assembly_width = static_cast<double>(shape_.first) * pitch_;
  after_baffle_ref_width_ = assembly_width - gap_width_ - baffle_width_;

  if (after_baffle_ref_width_ <= 0.) {
    auto mssg =
        "The gap width + baffle width is greater than the assembly width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

void PWRReflector::set_flux_tolerance(double ftol) {
  if (ftol <= 0.) {
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ftol >= 0.1) {
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  flux_tolerance_ = ftol;
}

void PWRReflector::set_keff_tolerance(double ktol) {
  if (ktol <= 0.) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ktol >= 0.1) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  keff_tolerance_ = ktol;
}

void PWRReflector::set_pins(const std::vector<Pin>& pins) {
  if (pins.size() != shape_.first * shape_.second) {
    auto mssg = "The number of pins does not agree with the assembly shape.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  pins_.reserve(pins.size());
  for (const auto& pin : pins) {
    if (std::holds_alternative<std::shared_ptr<FuelPin>>(pin)) {
      auto ptr = std::get<std::shared_ptr<FuelPin>>(pin);
      if (ptr == nullptr) {
        auto mssg = "All pins must be defined.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      pins_.push_back(std::make_shared<FuelPin>(*ptr));
    } else if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin)) {
      auto ptr = std::get<std::shared_ptr<GuideTube>>(pin);
      if (ptr == nullptr) {
        auto mssg = "All pins must be defined.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      pins_.push_back(std::make_shared<GuideTube>(*ptr));
    } else if (std::holds_alternative<std::shared_ptr<BurnablePoisonPin>>(
                   pin)) {
      auto ptr = std::get<std::shared_ptr<BurnablePoisonPin>>(pin);
      if (ptr == nullptr) {
        auto mssg = "All pins must be defined.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      pins_.push_back(std::make_shared<BurnablePoisonPin>(*ptr));
    } else {
      auto mssg = "Unsupported pin variant.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }
}

void PWRReflector::set_moderator(std::shared_ptr<Material> mod) {
  if (mod == nullptr) {
    auto mssg = "Moderator cannot be None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  moderator_ = mod;
  moderator_xs_ = moderator_->dilution_xs(
      std::vector<double>(moderator_->size(), 1.0E10), ndl_);
  moderator_xs_->set_name("Moderator");
}

void PWRReflector::set_num_azimuthal_angles(std::uint32_t n) {
  if (n < 4) {
    auto mssg = "Number of azimuthal angles must be >= 4.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (n % 2 != 0) {
    auto mssg = "Number of azimuthal angles must be even.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  num_azimuthal_angles_ = n;
}

void PWRReflector::set_track_spacing(double t) {
  if (t <= 0. || t >= 1.) {
    auto mssg = "Track spacing must be > 0 and < 1.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  track_spacing_ = t;
}

void PWRReflector::set_dancoff_track_spacing(double t) {
  if (t <= 0. || t >= 1.) {
    auto mssg = "Track spacing must be > 0 and < 1.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  dancoff_track_spacing_ = t;
}

void PWRReflector::set_dancoff_num_azimuthal_angles(std::uint32_t n) {
  if (n < 4) {
    auto mssg = "Number of azimuthal angles must be >= 4.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (n % 2 != 0) {
    auto mssg = "Number of azimuthal angles must be even.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  dancoff_num_azimuthal_angles_ = n;
}

void PWRReflector::solve() {
  if (pins_.size() == 0) {
    auto mssg = "Cannot solve PWR reflector problem. No pins provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (condensation_scheme_.size() == 0) {
    auto mssg =
        "Cannot solve PWR reflector problem. No energy condensation scheme "
        "provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (few_group_condensation_scheme_.size() == 0) {
    auto mssg =
        "Cannot solve PWR reflector problem. Few-group energy condensation "
        "scheme is empty.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  build_reflector_dancoff_geometry();
  get_fuel_dancoff_corrections();
  get_clad_dancoff_corrections();
  pin_cell_calc();
  condense_xs();
  baffle_spectrum_calc();
  moc_calc();
  few_group_xs();
  compute_adf_cdf();
}

double PWRReflector::isolated_fuel_pin_flux(DancoffMaterial dm) const {
  // We first make the system for an isolated fuel pin.
  // We isolate it by multiplying the pitch by 20.
  std::shared_ptr<SimplePinCell> isolated_fp{nullptr};
  for (const auto& pin : pins_) {
    if (std::holds_alternative<std::shared_ptr<FuelPin>>(pin)) {
      auto ptr = std::get<std::shared_ptr<FuelPin>>(pin);
      if (dm == DancoffMaterial::Fuel)
        isolated_fp = ptr->make_fuel_dancoff_cell(20. * pitch_, moderator_);
      else
        isolated_fp = ptr->make_clad_dancoff_cell(20. * pitch_, moderator_);
      break;
    }
  }

  if (isolated_fp == nullptr) {
    auto mssg = "No FuelPin type found in pins.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  std::shared_ptr<Cartesian2D> iso_geom = std::make_shared<Cartesian2D>(
      std::vector<double>{20. * pitch_}, std::vector<double>{20. * pitch_});
  iso_geom->set_tiles({isolated_fp});
  std::shared_ptr<MOCDriver> iso_moc = std::make_shared<MOCDriver>(iso_geom);

  // Set the source
  std::string mat_name;
  if (dm == DancoffMaterial::Fuel)
    mat_name = "Fuel";
  else
    mat_name = "Clad";
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() != mat_name) {
      // If we aren't in the fuel, set the source to be the value of
      // the potential xs (should be Et).
      iso_moc->set_extern_src(i, 0, i_xs->Et(0));
    }
  }

  // Solve the isolated pin problem
  iso_moc->x_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->x_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->generate_tracks(dancoff_num_azimuthal_angles_,
                           dancoff_track_spacing_, dancoff_polar_quadrature_);
  iso_moc->sim_mode() = SimulationMode::FixedSource;
  iso_moc->set_flux_tolerance(1.0e-5);
  iso_moc->solve();

  // Obtain isolated flux
  double iso_flux = 0.;
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() == mat_name) {
      iso_flux = iso_moc->flux(i, 0);
      break;
    }
  }

  return iso_flux;
}

double PWRReflector::isolated_guide_tube_flux() const {
  // We first make the system for an isolated guide tube.
  // We isolate it by multiplying the pitch by 20.
  std::shared_ptr<SimplePinCell> isolated_gt{nullptr};
  for (const auto& pin : pins_) {
    if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin)) {
      auto ptr = std::get<std::shared_ptr<GuideTube>>(pin);
      isolated_gt = ptr->make_clad_dancoff_cell(20. * pitch_, moderator_);
      break;
    }
  }

  if (isolated_gt == nullptr) {
    return 0.;
  }

  std::shared_ptr<Cartesian2D> iso_geom = std::make_shared<Cartesian2D>(
      std::vector<double>{20. * pitch_}, std::vector<double>{20. * pitch_});
  iso_geom->set_tiles({isolated_gt});
  std::shared_ptr<MOCDriver> iso_moc = std::make_shared<MOCDriver>(iso_geom);

  // Set the source
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() != "Clad") {
      // If we aren't in the fuel, set the source to be the value of
      // the potential xs (should be Et).
      iso_moc->set_extern_src(i, 0, i_xs->Et(0));
    }
  }

  // Solve the isolated pin problem
  iso_moc->x_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->x_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->generate_tracks(dancoff_num_azimuthal_angles_,
                           dancoff_track_spacing_, dancoff_polar_quadrature_);
  iso_moc->sim_mode() = SimulationMode::FixedSource;
  iso_moc->set_flux_tolerance(1.0e-5);
  iso_moc->solve();

  // Obtain isolated flux
  double iso_flux = 0.;
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() == "Clad") {
      iso_flux = iso_moc->flux(i, 0);
      break;
    }
  }

  return iso_flux;
}

double PWRReflector::isolated_burnable_poison_tube_flux() const {
  // We first make the system for an isolated burnable poison pin.
  // We isolate it by multiplying the pitch by 20.
  std::shared_ptr<SimplePinCell> isolated_bp{nullptr};
  for (const auto& pin : pins_) {
    if (std::holds_alternative<std::shared_ptr<BurnablePoisonPin>>(pin)) {
      auto ptr = std::get<std::shared_ptr<BurnablePoisonPin>>(pin);
      isolated_bp = ptr->make_clad_dancoff_cell(20. * pitch_, moderator_);
      break;
    }
  }

  if (isolated_bp == nullptr) {
    return 0.;
  }

  std::shared_ptr<Cartesian2D> iso_geom = std::make_shared<Cartesian2D>(
      std::vector<double>{20. * pitch_}, std::vector<double>{20. * pitch_});
  iso_geom->set_tiles({isolated_bp});
  std::shared_ptr<MOCDriver> iso_moc = std::make_shared<MOCDriver>(iso_geom);

  // Set the source
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() != "Clad") {
      // If we aren't in the fuel, set the source to be the value of
      // the potential xs (should be Et).
      iso_moc->set_extern_src(i, 0, i_xs->Et(0));
    }
  }

  // Solve the isolated pin problem
  iso_moc->x_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->x_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_min_bc() = BoundaryCondition::Vacuum;
  iso_moc->y_max_bc() = BoundaryCondition::Vacuum;
  iso_moc->generate_tracks(dancoff_num_azimuthal_angles_,
                           dancoff_track_spacing_, dancoff_polar_quadrature_);
  iso_moc->sim_mode() = SimulationMode::FixedSource;
  iso_moc->set_flux_tolerance(1.0e-5);
  iso_moc->solve();

  // Obtain isolated flux
  double iso_flux = 0.;
  for (std::size_t i = 0; i < iso_moc->nfsr(); i++) {
    const auto& i_xs = iso_moc->xs(i);
    if (i_xs->name() == "Clad") {
      iso_flux = iso_moc->flux(i, 0);
      break;
    }
  }

  return iso_flux;
}

void PWRReflector::build_reflector_dancoff_geometry() {
  // We start by making cross sections for the different materials.
  xt::xtensor<double, 1> Et = {moderator_->potential_xs()};
  xt::xtensor<double, 1> Ea = {Et(0)};
  xt::xtensor<double, 2> Es = {{0.}};
  auto mod = std::make_shared<CrossSection>(Et, Ea, Es, "Moderator");

  std::shared_ptr<CrossSection> baff{nullptr};
  if (baffle_) {
    Et(0) = baffle_->potential_xs();
    Ea(0) = Et(0);
    baff = std::make_shared<CrossSection>(Et, Ea, Es, "Baffle");
  }

  // Determine number of regions
  std::size_t NG = static_cast<std::size_t>(gap_width_ / 0.3);
  if (NG == 0) NG++;
  double dx_gap = gap_width_ / static_cast<double>(NG);

  std::size_t NB = 0;
  double dx_baffle = 0.;
  if (baffle_) {
    NB = static_cast<std::size_t>(baffle_width_ / 0.3);
    dx_baffle = baffle_width_ / static_cast<double>(NB);
  }

  std::size_t NR = static_cast<std::size_t>(after_baffle_ref_width_ / 0.3) + 1;
  double dx_ref = after_baffle_ref_width_ / static_cast<double>(NR);

  const double asmbly_y = static_cast<double>(shape_.second) * pitch_;
  std::size_t NY = static_cast<std::size_t>(asmbly_y / 0.3) + 1;
  const double delta_y = asmbly_y / static_cast<double>(NY);

  // Create base tiles
  std::shared_ptr<EmptyCell> gap_tile =
      std::make_shared<EmptyCell>(mod, dx_gap, delta_y);
  std::shared_ptr<EmptyCell> baff_tile{nullptr};
  if (baffle_)
    baff_tile = std::make_shared<EmptyCell>(baff, dx_baffle, delta_y);
  std::shared_ptr<EmptyCell> ref_tile =
      std::make_shared<EmptyCell>(mod, dx_ref, delta_y);

  // Create the dx and dy arrays
  std::vector<double> dy(NY, delta_y);
  std::vector<double> dx;
  dx.reserve(NG + NB + NR);
  for (std::size_t i = 0; i < NG + NB + NR; i++) {
    if (i < NG)
      dx.push_back(dx_gap);
    else if (i < NG + NB)
      dx.push_back(dx_baffle);
    else
      dx.push_back(dx_ref);
  }

  reflector_dancoff_geom_ = std::make_shared<Cartesian2D>(dx, dy);

  std::vector<Cartesian2D::TileFill> tiles;
  tiles.reserve((NG + NB + NR) * NY);
  for (std::size_t j = 0; j < NY; j++) {
    for (std::size_t i = 0; i < NG; i++) tiles.push_back(gap_tile);
    for (std::size_t i = 0; i < NB; i++) tiles.push_back(baff_tile);
    for (std::size_t i = 0; i < NR; i++) tiles.push_back(ref_tile);
  }

  reflector_dancoff_geom_->set_tiles(tiles);
}

void PWRReflector::build_reflector_geometry() {
  // Determine number of regions
  std::size_t NG = static_cast<std::size_t>(gap_width_ / 0.3);
  if (NG == 0) NG++;
  double dx_gap = gap_width_ / static_cast<double>(NG);

  std::size_t NB = 0;
  double dx_baffle = 0.;
  if (baffle_) {
    NB = static_cast<std::size_t>(baffle_width_ / 0.3);
    dx_baffle = baffle_width_ / static_cast<double>(NB);
  }

  std::size_t NR = static_cast<std::size_t>(after_baffle_ref_width_ / 0.3) + 1;
  double dx_ref = after_baffle_ref_width_ / static_cast<double>(NR);

  const double asmbly_y = static_cast<double>(shape_.second) * pitch_;
  std::size_t NY = static_cast<std::size_t>(asmbly_y / 0.3) + 1;
  const double delta_y = asmbly_y / static_cast<double>(NY);

  // Create base tiles
  std::shared_ptr<EmptyCell> gap_tile =
      std::make_shared<EmptyCell>(macro_gap_xs_, dx_gap, delta_y);
  std::shared_ptr<EmptyCell> baff_tile{nullptr};
  if (baffle_)
    baff_tile =
        std::make_shared<EmptyCell>(macro_baffle_xs_, dx_baffle, delta_y);
  std::shared_ptr<EmptyCell> ref_tile =
      std::make_shared<EmptyCell>(macro_ref_xs_, dx_ref, delta_y);

  // Create the dx and dy arrays
  std::vector<double> dy(NY, delta_y);
  std::vector<double> dx;
  dx.reserve(NG + NB + NR);
  for (std::size_t i = 0; i < NG + NB + NR; i++) {
    if (i < NG)
      dx.push_back(dx_gap);
    else if (i < NG + NB)
      dx.push_back(dx_baffle);
    else
      dx.push_back(dx_ref);
  }

  moc_refl_geom_ = std::make_shared<Cartesian2D>(dx, dy);

  std::vector<Cartesian2D::TileFill> tiles;
  tiles.reserve((NG + NB + NR) * NY);
  for (std::size_t j = 0; j < NY; j++) {
    for (std::size_t i = 0; i < NG; i++) tiles.push_back(gap_tile);
    for (std::size_t i = 0; i < NB; i++) tiles.push_back(baff_tile);
    for (std::size_t i = 0; i < NR; i++) tiles.push_back(ref_tile);
  }

  moc_refl_geom_->set_tiles(tiles);
}

void PWRReflector::get_fuel_dancoff_corrections() {
  spdlog::info("");
  spdlog::info("Computing Dancoff factors for fuel");
  set_logging_level(LogLevel::warn);

  // First get the flux of the isolated pin
  const double iso_flux = isolated_fuel_pin_flux(DancoffMaterial::Fuel);

  // Now we setup the lattice problem
  std::vector<Cartesian2D::TileFill> fuel_df_pins;
  fuel_df_pins.reserve(pins_.size());
  for (const auto& pin : pins_) {
    fuel_df_pins.push_back(std::visit(
        [this](const auto& P) {
          return P->make_fuel_dancoff_cell(pitch_, moderator_);
        },
        pin));
  }

  std::shared_ptr<Cartesian2D> asmbly_geom =
      std::make_shared<Cartesian2D>(std::vector<double>(shape_.first, pitch_),
                                    std::vector<double>(shape_.first, pitch_));
  asmbly_geom->set_tiles(fuel_df_pins);

  const double asmbly_dx = static_cast<double>(shape_.first) * pitch_;
  const double asmbly_dy = static_cast<double>(shape_.second) * pitch_;
  std::shared_ptr<Cartesian2D> geom = std::make_shared<Cartesian2D>(
      std::vector<double>(2, asmbly_dx), std::vector<double>(1, asmbly_dy));
  std::vector<Cartesian2D::TileFill> tiles{asmbly_geom,
                                           reflector_dancoff_geom_};
  geom->set_tiles(tiles);
  std::shared_ptr<MOCDriver> moc = std::make_shared<MOCDriver>(geom);

  // Set the source
  for (std::size_t i = 0; i < moc->nfsr(); i++) {
    const auto& i_xs = moc->xs(i);
    if (i_xs->name() != "Fuel") {
      // If we aren't in the fuel, set the source to be the value of
      // the potential xs (should be Et).
      moc->set_extern_src(i, 0, i_xs->Et(0));
    }
  }

  // Solve the lattice problem
  moc->x_max_bc() = BoundaryCondition::Vacuum;
  moc->generate_tracks(dancoff_num_azimuthal_angles_, dancoff_track_spacing_,
                       dancoff_polar_quadrature_);
  moc->sim_mode() = SimulationMode::FixedSource;
  moc->set_flux_tolerance(1.0e-5);
  moc->solve();

  // Now we need to calculate the dancoff correction for each pin
  fuel_dancoff_corrections_.reserve(pins_.size());
  const Direction u(1.0, 0.0);
  for (std::size_t j = 0; j < shape_.second; j++) {
    const double y = moc->y_max() - (static_cast<double>(j) + 0.5) * pitch_;
    for (std::size_t i = 0; i < shape_.first; i++) {
      const double x = moc->x_min() + (static_cast<double>(i) + 0.5) * pitch_;
      const Vector r(x, y);
      const auto& xs = moc->xs(r, u);
      if (xs->name() == "Fuel") {
        const double flux = moc->flux(r, u, 0);
        const double C = (iso_flux - flux) / iso_flux;
        fuel_dancoff_corrections_.push_back(C);
      } else {
        fuel_dancoff_corrections_.push_back(0.0);
      }
    }
  }

  set_logging_level(LogLevel::info);
}

void PWRReflector::get_clad_dancoff_corrections() {
  spdlog::info("");
  spdlog::info("Computing Dancoff factors for cladding");
  set_logging_level(LogLevel::warn);

  // obtain isolated fluxes
  const double iso_flux_fp = isolated_fuel_pin_flux(DancoffMaterial::Clad);
  const double iso_flux_gt = isolated_guide_tube_flux();
  const double iso_flux_bp = isolated_burnable_poison_tube_flux();

  // Now we setup the lattice problem
  std::vector<Cartesian2D::TileFill> fuel_df_pins;
  fuel_df_pins.reserve(pins_.size());
  for (const auto& pin : pins_) {
    fuel_df_pins.push_back(std::visit(
        [this](const auto& P) {
          return P->make_clad_dancoff_cell(pitch_, moderator_);
        },
        pin));
  }

  std::shared_ptr<Cartesian2D> asmbly_geom =
      std::make_shared<Cartesian2D>(std::vector<double>(shape_.first, pitch_),
                                    std::vector<double>(shape_.first, pitch_));
  asmbly_geom->set_tiles(fuel_df_pins);

  const double asmbly_dx = static_cast<double>(shape_.first) * pitch_;
  const double asmbly_dy = static_cast<double>(shape_.second) * pitch_;
  std::shared_ptr<Cartesian2D> geom = std::make_shared<Cartesian2D>(
      std::vector<double>(2, asmbly_dx), std::vector<double>(1, asmbly_dy));
  std::vector<Cartesian2D::TileFill> tiles{asmbly_geom,
                                           reflector_dancoff_geom_};
  geom->set_tiles(tiles);
  std::shared_ptr<MOCDriver> moc = std::make_shared<MOCDriver>(geom);

  // Set the source
  for (std::size_t i = 0; i < moc->nfsr(); i++) {
    const auto& i_xs = moc->xs(i);
    if (i_xs->name() != "Clad") {
      // If we aren't in the fuel, set the source to be the value of
      // the potential xs (should be Et).
      moc->set_extern_src(i, 0, i_xs->Et(0));
    }
  }

  // Solve the lattice problem
  moc->x_max_bc() = BoundaryCondition::Vacuum;
  moc->generate_tracks(dancoff_num_azimuthal_angles_, dancoff_track_spacing_,
                       dancoff_polar_quadrature_);
  moc->sim_mode() = SimulationMode::FixedSource;
  moc->set_flux_tolerance(1.0e-5);
  moc->solve();

  // Now we need to calculate the dancoff correction for each pin
  clad_dancoff_corrections_.reserve(pins_.size());
  const Direction u(1.0, 0.0);
  std::size_t i_pin = 0;
  for (std::size_t j = 0; j < shape_.second; j++) {
    const double y = moc->y_max() - (static_cast<double>(j) + 0.5) * pitch_;
    for (std::size_t i = 0; i < shape_.first; i++) {
      const double x = moc->x_min() + (static_cast<double>(i) + 0.5) * pitch_;
      const auto& pin = pins_[i_pin];
      const Vector r =
          Vector(x, y) +
          std::visit([](const auto& P) { return P->clad_offset(); }, pin);
      const auto& xs = moc->xs(r, u);
      const double flux = moc->flux(r, u, 0);

      double C = 0.;
      if (xs->name() == "Clad") {
        if (std::holds_alternative<std::shared_ptr<FuelPin>>(pin)) {
          C = (iso_flux_fp - flux) / iso_flux_fp;
        } else if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin)) {
          C = (iso_flux_gt - flux) / iso_flux_gt;
        } else if (std::holds_alternative<std::shared_ptr<BurnablePoisonPin>>(
                       pin)) {
          C = (iso_flux_bp - flux) / iso_flux_bp;
        }
      }
      clad_dancoff_corrections_.push_back(C);
      i_pin++;
    }
  }

  set_logging_level(LogLevel::info);
}

void PWRReflector::pin_cell_calc() {
  spdlog::info("");
  spdlog::info("Performing micro-group pin cell calcuations");
  spdlog::info("Please wait...");
  set_logging_level(LogLevel::warn);

  // Preload all nuclides
  moderator_->load_nuclides(ndl_);
  for (const auto& pin : pins_)
    std::visit([this](const auto& P) { P->load_nuclides(ndl_); }, pin);

  // First, get all fuel pin indices
  std::vector<std::size_t> fp_inds;
  std::vector<std::size_t> other_inds;
  fp_inds.reserve(pins_.size());
  other_inds.reserve(pins_.size());
  for (std::size_t i = 0; i < pins_.size(); i++) {
    if (std::holds_alternative<std::shared_ptr<FuelPin>>(pins_[i])) {
      fp_inds.push_back(i);
    } else {
      other_inds.push_back(i);
    }
  }

  // Resize the cell and flux arrays
  pin_1d_cells.resize(pins_.size(), nullptr);
  pin_1d_fluxes.resize(pins_.size(), nullptr);

  // Now do all fuel pins
#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(fp_inds.size()); ii++) {
    const std::size_t i = fp_inds[static_cast<std::size_t>(ii)];
    const auto& pin = std::get<std::shared_ptr<FuelPin>>(pins_[i]);
    const double fuel_dancoff = fuel_dancoff_corrections_[i];
    const double clad_dancoff = clad_dancoff_corrections_[i];

    pin_1d_cells[i] = pin->make_cylindrical_cell(
        pitch_, fuel_dancoff, moderator_xs_, ndl_, clad_dancoff);
    pin_1d_cells[i]->solve();
    pin_1d_fluxes[i] = std::make_shared<CylindricalFluxSolver>(pin_1d_cells[i]);
    pin_1d_fluxes[i]->solve();
  }

  // Compute the homogenized xs for the average fuel pin cell
  avg_fp_ = pin_1d_fluxes[fp_inds[0]]->homogenize();
  for (std::size_t i = 1; i < fp_inds.size(); i++) {
    *avg_fp_ += *(pin_1d_fluxes[fp_inds[i]]->homogenize());
  }
  *avg_fp_ *= 1. / static_cast<double>(fp_inds.size());

  // Now we compute all of the non-fuel cells
  const double buffer_rad = std::sqrt(9. * pitch_ * pitch_ / PI);
#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(other_inds.size()); ii++) {
    const std::size_t i = other_inds[static_cast<std::size_t>(ii)];
    const double clad_dancoff = clad_dancoff_corrections_[i];

    const auto& pin_var = pins_[i];
    if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin_var)) {
      const auto& pin = std::get<std::shared_ptr<GuideTube>>(pin_var);
      pin_1d_cells[i] = pin->make_cylindrical_cell(
          pitch_, moderator_xs_, buffer_rad, avg_fp_, ndl_, clad_dancoff);
    } else if (std::holds_alternative<std::shared_ptr<BurnablePoisonPin>>(
                   pin_var)) {
      const auto& pin = std::get<std::shared_ptr<BurnablePoisonPin>>(pin_var);
      pin_1d_cells[i] = pin->make_cylindrical_cell(
          pitch_, moderator_xs_, buffer_rad, avg_fp_, ndl_, clad_dancoff);
    }

    pin_1d_cells[i]->solve();
    pin_1d_fluxes[i] = std::make_shared<CylindricalFluxSolver>(pin_1d_cells[i]);
    pin_1d_fluxes[i]->solve();
  }

  // At this point, we no longer need the raw nuclear data, so we unload it
  // to conserve memory
  ndl_->unload();

  set_logging_level(LogLevel::info);
}

void PWRReflector::baffle_spectrum_calc() {
  spdlog::info("");
  spdlog::info("Performing the micro-group baffle/reflector calcuation");
  spdlog::info("Please wait...");
  set_logging_level(LogLevel::warn);

  // First, we create the baffle cross section (if present)
  if (baffle_) {
    const double Ee = 1. / (2. * baffle_width_);
    baffle_xs_ = baffle_->roman_xs(0., Ee, ndl_);
    baffle_xs_->set_name("Baffle");
  }

  // Now we create a 1D annular pin cell problem to approximate the fine
  // specrum in the gap, baffle, and reflector
  std::size_t NF = 2 * 17;
  std::vector<std::size_t> fuel_regions(NF);
  std::iota(fuel_regions.begin(), fuel_regions.end(), 0);
  double dx_fuel = pitch_;

  std::size_t NG = static_cast<std::size_t>(gap_width_ / 0.3);
  if (NG == 0) NG++;
  double dx_gap = gap_width_ / static_cast<double>(NG);
  std::vector<std::size_t> gap_regions(NG);
  std::iota(gap_regions.begin(), gap_regions.end(), NF);

  std::size_t NB = 0;
  double dx_baffle = 0.;
  std::vector<std::size_t> baffle_regions;
  if (baffle_) {
    NB = static_cast<std::size_t>(baffle_width_ / 0.3);
    dx_baffle = baffle_width_ / static_cast<double>(NB);
    baffle_regions.resize(NB);
    std::iota(baffle_regions.begin(), baffle_regions.end(), NF + NG);
  }

  std::size_t NR = static_cast<std::size_t>(after_baffle_ref_width_ / 0.3) + 1;
  double dx_ref = after_baffle_ref_width_ / static_cast<double>(NR);
  std::vector<std::size_t> ref_regions(NR);
  std::iota(ref_regions.begin(), ref_regions.end(), NF + NG + NB);

  std::vector<double> radii;
  std::vector<std::shared_ptr<CrossSection>> mats;
  radii.reserve(NF + NG + NB + NR);
  mats.reserve(NF + NG + NB + NR);

  radii.push_back(dx_fuel);
  mats.push_back(avg_fp_);
  for (std::size_t i = 1; i < NF; i++) {
    radii.push_back(radii.back() + dx_fuel);
    mats.push_back(avg_fp_);
  }
  for (std::size_t i = 0; i < NG; i++) {
    radii.push_back(radii.back() + dx_gap);
    mats.push_back(moderator_xs_);
  }
  for (std::size_t i = 0; i < NB; i++) {
    radii.push_back(radii.back() + dx_baffle);
    mats.push_back(baffle_xs_);
  }
  for (std::size_t i = 0; i < NR; i++) {
    radii.push_back(radii.back() + dx_ref);
    mats.push_back(moderator_xs_);
  }

  reflector_cyl_cell_ = std::make_shared<CylindricalCell>(radii, mats);
  reflector_cyl_cell_->solve(true); // Solve in parallel
  reflector_cyl_flux_cell_ =
      std::make_shared<CylindricalFluxSolver>(reflector_cyl_cell_);
  reflector_cyl_flux_cell_->set_albedo(0.);
  reflector_cyl_flux_cell_->solve(true); // Solve in parallel

  // Now we condense the macro cross sections
  auto gap_homog_xs = reflector_cyl_flux_cell_->homogenize(gap_regions);
  auto gap_spectrum =
      reflector_cyl_flux_cell_->homogenize_flux_spectrum(gap_regions);
  macro_gap_xs_ = gap_homog_xs->condense(condensation_scheme_, gap_spectrum);

  auto ref_homog_xs = reflector_cyl_flux_cell_->homogenize(ref_regions);
  auto ref_spectrum =
      reflector_cyl_flux_cell_->homogenize_flux_spectrum(ref_regions);
  macro_ref_xs_ = ref_homog_xs->condense(condensation_scheme_, ref_spectrum);

  if (baffle_) {
    auto baf_homog_xs = reflector_cyl_flux_cell_->homogenize(baffle_regions);
    auto baf_spectrum =
        reflector_cyl_flux_cell_->homogenize_flux_spectrum(baffle_regions);
    macro_baffle_xs_ =
        baf_homog_xs->condense(condensation_scheme_, baf_spectrum);
  }

  set_logging_level(LogLevel::info);
}

void PWRReflector::condense_xs() {
  spdlog::info("");
  spdlog::info("Performing pin cell energy condensation");
  set_logging_level(LogLevel::warn);

  for (int ii = 0; ii < static_cast<int>(pins_.size()); ii++) {
    const std::size_t i = static_cast<std::size_t>(ii);
    auto& pin = pins_[i];
    const auto& cell_flux = pin_1d_fluxes[i];

    std::size_t NR = cell_flux->nregions();
    if (std::holds_alternative<std::shared_ptr<FuelPin>>(pin) == false) NR--;

    for (std::size_t r = 0; r < NR; r++) {
      const auto& xs = cell_flux->xs(r);
      auto flux_spectrum =
          cell_flux->homogenize_flux_spectrum(std::vector<std::size_t>{r});
      std::visit(
          [this, &xs, &flux_spectrum](auto& P) {
            P->condensed_xs().push_back(
                xs->condense(condensation_scheme_, flux_spectrum));
            P->condensed_xs().back()->set_name(xs->name());
          },
          pin);
    }
  }

  set_logging_level(LogLevel::info);
}

void PWRReflector::moc_calc() {
  spdlog::info("");
  spdlog::info("Performing macrogroup assembly calculation");

  std::vector<Cartesian2D::TileFill> moc_pins;
  moc_pins.reserve(pins_.size());
  for (std::size_t i = 0; i < pins_.size(); i++) {
    const auto& pin = pins_[i];
    moc_pins.push_back(std::visit(
        [this](const auto& P) { return P->make_moc_cell(pitch_); }, pin));
  }

  std::vector<double> dx(shape_.first, pitch_);
  std::vector<double> dy(shape_.second, pitch_);
  moc_asmbly_geom_ = std::make_shared<Cartesian2D>(dx, dy);
  moc_asmbly_geom_->set_tiles(moc_pins);

  build_reflector_geometry();

  const double asmbly_dx = static_cast<double>(shape_.first) * pitch_;
  const double asmbly_dy = static_cast<double>(shape_.second) * pitch_;
  moc_geom_ = std::make_shared<Cartesian2D>(std::vector<double>(2, asmbly_dx),
                                            std::vector<double>(1, asmbly_dy));
  std::vector<Cartesian2D::TileFill> tiles{moc_asmbly_geom_, moc_refl_geom_};
  moc_geom_->set_tiles(tiles);

  moc_ = std::make_shared<MOCDriver>(moc_geom_);
  if (plot_assembly_) {
    ImApp::App guiplotter(1920, 1080, "Scarabee MOC Plotter");
    guiplotter.enable_docking();
    guiplotter.push_layer(std::make_unique<MOCPlotter>(moc_.get()));
    guiplotter.run();
  }

  moc_->x_max_bc() = BoundaryCondition::Vacuum;
  moc_->generate_tracks(num_azimuthal_angles_, track_spacing_,
                        polar_quadrature_);
  moc_->set_keff_tolerance(keff_tolerance_);
  moc_->set_flux_tolerance(flux_tolerance_);
  moc_->solve();
}

void PWRReflector::few_group_xs() {
  spdlog::info("");
  spdlog::info("Generating few group cross sections");

  std::vector<std::size_t> asmbly_regions(moc_asmbly_geom_->num_fsrs());
  std::iota(asmbly_regions.begin(), asmbly_regions.end(), 0);
  std::vector<std::size_t> refl_regions(moc_refl_geom_->num_fsrs());
  std::iota(refl_regions.begin(), refl_regions.end(), asmbly_regions.size());

  asmbly_diffusion_xs_ = make_diffusion_xs(asmbly_regions);
  refl_diffusion_xs_ = make_diffusion_xs(refl_regions);

  auto NG = refl_diffusion_xs_->ngroups();
  xt::xtensor<double, 1> D = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 2> Es = xt::zeros<double>({NG, NG});

  for (std::size_t g = 0; g < NG; g++) {
    D(g) = refl_diffusion_xs_->D(g);
    Ea(g) = refl_diffusion_xs_->Ea(g);
    for (std::size_t gg = 0; gg < NG; gg++) {
      Es(g, gg) = refl_diffusion_xs_->Es(g, gg);
    }
  }

  std::stringstream D_str, Ea_str, Es_str;
  D_str << "D  : " << D;
  Ea_str << "Ea : " << Ea;
  spdlog::info(D_str.str());
  spdlog::info(Ea_str.str());

  Es_str << "Es : " << xt::view(Es, 0, xt::all());
  spdlog::info(Es_str.str());
  for (std::size_t g = 1; g < NG; g++) {
    Es_str.str(std::string());
    Es_str << "     " << xt::view(Es, g, xt::all());
    spdlog::info(Es_str.str());
  }
}

std::shared_ptr<DiffusionCrossSection> PWRReflector::make_diffusion_xs(
    const std::vector<std::size_t>& regions) const {
  const auto homog_xs = moc_->homogenize(regions);
  const auto flux_spectrum = moc_->homogenize_flux_spectrum(regions);
  auto NG = homog_xs->ngroups();
  const bool fissile = homog_xs->fissile();

  xt::xtensor<double, 1> D = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 2> Es = xt::zeros<double>({NG, NG});
  xt::xtensor<double, 1> Ef, vEf, chi;
  if (fissile) {
    Ef = xt::zeros<double>({NG});
    vEf = xt::zeros<double>({NG});
    chi = xt::zeros<double>({NG});
  }

  for (std::size_t g = 0; g < NG; g++) {
    D(g) = 1. / (3. * homog_xs->Etr(g));
    Ea(g) = homog_xs->Ea(g);

    if (fissile) {
      Ef(g) = homog_xs->Ef(g);
      vEf(g) = homog_xs->vEf(g);
      chi(g) = homog_xs->chi(g);
    }

    for (std::size_t gg = 0; gg < NG; gg++) {
      Es(g, gg) = homog_xs->Es_tr(g, gg);
    }
  }

  std::unique_ptr<DiffusionCrossSection> diff_xs{nullptr};
  if (fissile) {
    diff_xs = std::make_unique<DiffusionCrossSection>(D, Ea, Es, Ef, vEf, chi);
  } else {
    diff_xs = std::make_unique<DiffusionCrossSection>(D, Ea, Es);
  }

  // According to Smith, one should do energy condensation on the
  // diffusion coefficients, and not on the transport cross sections which
  // one could then use to make diffusion coefficients [1]. This is in
  // contradiction to Lattice Physics Computations which states that
  // either method is acceptable [2]. In light of these comments, I have
  // chosen to go with Smith's recommendation of performing energy
  // condensation on the diffusion coefficients.
  return diff_xs->condense(few_group_condensation_scheme_, flux_spectrum);
}

void PWRReflector::compute_adf_cdf() {}

void PWRReflector::save_diffusion_data(const std::string& fname) const {
  if (refl_diffusion_xs_ == nullptr || adf_.size() == 0 || cdf_.size() == 0) {
    auto mssg = "Cannot save DiffusionData. Assembly has not been solved.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  DiffusionData dd(refl_diffusion_xs_);
  dd.set_adf(adf_);
  dd.set_cdf(cdf_);
  dd.save(fname);
}

std::vector<double> PWRReflector::compute_avg_surface_flx(
    const std::vector<std::pair<std::size_t, double>>& segments) const {
  double total_len = 0.;
  for (const auto& i_l : segments) total_len += i_l.second;
  const double invs_tot_len = 1. / total_len;

  // Number of few groups
  const std::size_t NG = few_group_condensation_scheme_.size();

  std::vector<double> flx(NG, 0.);

  for (const auto& i_l : segments) {
    for (std::size_t G = 0; G < NG; G++) {
      const std::size_t gmin = few_group_condensation_scheme_[G].first;
      const std::size_t gmax = few_group_condensation_scheme_[G].first;

      for (std::size_t g = gmin; g <= gmax; g++) {
        flx[G] += i_l.second * moc_->flux(i_l.first, g);
      }
    }
  }

  for (auto& f : flx) f *= invs_tot_len;

  return flx;
}

std::vector<double> PWRReflector::compute_avg_flx(const Vector& r,
                                                  const Direction& u) const {
  const auto fsr = moc_->get_fsr(r, u);
  const std::size_t i = moc_->get_fsr_indx(fsr);

  // Number of few groups
  const std::size_t NG = few_group_condensation_scheme_.size();

  std::vector<double> flx(NG, 0.);

  for (std::size_t G = 0; G < NG; G++) {
    const std::size_t gmin = few_group_condensation_scheme_[G].first;
    const std::size_t gmax = few_group_condensation_scheme_[G].first;

    for (std::size_t g = gmin; g <= gmax; g++) {
      flx[G] += moc_->flux(i, g);
    }
  }

  return flx;
}

}  // namespace scarabee

// REFERENCES
// [1] K. S. Smith, “Nodal diffusion methods and lattice physics data in LWR
//     analyses: Understanding numerous subtle details,” Prog Nucl Energ,
//     vol. 101, pp. 360–369, 2017, doi: 10.1016/j.pnucene.2017.06.013.
// [2] D. Knott and A. Yamamoto, "Lattice Physics Computations" in
//     Handbook of Nuclear Engineering, 2010, p 1226.
