#include <assemblies/pwr_assembly.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/criticality_spectrum.hpp>
#include <utils/timer.hpp>
#include <moc/moc_plotter.hpp>
#include <diffusion/diffusion_data.hpp>

#include <xtensor/io/xio.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

namespace scarabee {

PWRAssembly::PWRAssembly(double pitch, std::shared_ptr<Material> moderator,
                         std::pair<std::size_t, std::size_t> shape,
                         std::shared_ptr<NDLibrary> ndl)
    : pitch_(pitch), shape_(shape), ndl_(ndl) {
  this->set_moderator(moderator);

  if (pitch_ <= 0.) {
    auto mssg = "Pitch must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ndl_->macro_group_condensation_scheme()) {
    this->set_condensation_scheme(
        ndl_->macro_group_condensation_scheme().value());
  }

  if (ndl_->few_group_condensation_scheme()) {
    this->set_few_group_condensation_scheme(
        ndl_->few_group_condensation_scheme().value());
  }
}

void PWRAssembly::set_flux_tolerance(double ftol) {
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

void PWRAssembly::set_keff_tolerance(double ktol) {
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

void PWRAssembly::set_criticality_spectrum_method(
    const std::optional<std::string>& csm) {
  if (csm.has_value() != false && *csm != "B1" && *csm != "b1" &&
      *csm != "P1" && *csm != "p1") {
    auto mssg = "Unknown criticality spectrum method.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (csm.has_value() == false) {
    criticality_spectrum_method_ = csm;
  } else if (*csm == "B1" || *csm == "b1") {
    criticality_spectrum_method_ = "B1";
  } else {
    criticality_spectrum_method_ = "P1";
  }
}

void PWRAssembly::set_pins(const std::vector<Pin>& pins) {
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
      pins_.push_back(ptr->clone());
    } else if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin)) {
      auto ptr = std::get<std::shared_ptr<GuideTube>>(pin);
      if (ptr == nullptr) {
        auto mssg = "All pins must be defined.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      pins_.push_back(ptr->clone());
    } else if (std::holds_alternative<std::shared_ptr<BurnablePoisonPin>>(
                   pin)) {
      auto ptr = std::get<std::shared_ptr<BurnablePoisonPin>>(pin);
      if (ptr == nullptr) {
        auto mssg = "All pins must be defined.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      pins_.push_back(ptr->clone());
    } else {
      auto mssg = "Unsupported pin variant.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }
}

void PWRAssembly::set_moderator(std::shared_ptr<Material> mod) {
  if (mod == nullptr) {
    auto mssg = "Moderator cannot be None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  moderator_ = mod;
  moderator_xs_ = moderator_->dilution_xs(
      std::vector<double>(moderator_->size(), 1.0E10), ndl_);

  if (moderator_xs_->name() == "") moderator_xs_->set_name("Moderator");
}

void PWRAssembly::set_num_azimuthal_angles(std::uint32_t n) {
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

void PWRAssembly::set_track_spacing(double t) {
  if (t <= 0. || t >= 1.) {
    auto mssg = "Track spacing must be > 0 and < 1.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  track_spacing_ = t;
}

void PWRAssembly::set_dancoff_track_spacing(double t) {
  if (t <= 0. || t >= 1.) {
    auto mssg = "Track spacing must be > 0 and < 1.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  dancoff_track_spacing_ = t;
}

void PWRAssembly::set_dancoff_num_azimuthal_angles(std::uint32_t n) {
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

void PWRAssembly::solve() {
  if (ndl_ == nullptr) {
    auto mssg =
        "Cannot solve PWR assembly problem. Nuclear data library is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (pins_.size() == 0) {
    auto mssg = "Cannot solve PWR assembly problem. No pins provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (condensation_scheme_.size() == 0) {
    auto mssg =
        "Cannot solve PWR assembly problem. No energy condensation scheme "
        "provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (few_group_condensation_scheme_.size() == 0) {
    auto mssg =
        "Cannot solve PWR assembly problem. Few-group energy condensation "
        "scheme is empty.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  get_fuel_dancoff_corrections();
  get_clad_dancoff_corrections();
  pin_cell_calc();
  condense_xs();
  moc_calc();
  criticality_spectrum();
  compute_form_factors();
  few_group_xs();
  compute_adf_cdf();
}

double PWRAssembly::isolated_fuel_pin_flux(DancoffMaterial dm) const {
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

double PWRAssembly::isolated_guide_tube_flux() const {
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

double PWRAssembly::isolated_burnable_poison_tube_flux() const {
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

void PWRAssembly::get_fuel_dancoff_corrections() {
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
  std::shared_ptr<Cartesian2D> geom =
      std::make_shared<Cartesian2D>(std::vector<double>(shape_.first, pitch_),
                                    std::vector<double>(shape_.first, pitch_));
  geom->set_tiles(fuel_df_pins);
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
  moc->x_min_bc() = this->boundary_conditions();
  moc->x_max_bc() = this->boundary_conditions();
  moc->y_min_bc() = this->boundary_conditions();
  moc->y_max_bc() = this->boundary_conditions();
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
        fuel_dancoff_corrections_.push_back(0.);
      }
    }
  }

  set_logging_level(LogLevel::info);
}

void PWRAssembly::get_clad_dancoff_corrections() {
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
  std::shared_ptr<Cartesian2D> geom =
      std::make_shared<Cartesian2D>(std::vector<double>(shape_.first, pitch_),
                                    std::vector<double>(shape_.first, pitch_));
  geom->set_tiles(fuel_df_pins);
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
  moc->x_min_bc() = this->boundary_conditions();
  moc->x_max_bc() = this->boundary_conditions();
  moc->y_min_bc() = this->boundary_conditions();
  moc->y_max_bc() = this->boundary_conditions();
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

void PWRAssembly::pin_cell_calc() {
  Timer pin_cell_calc_timer;
  pin_cell_calc_timer.start();

  spdlog::info("");
  spdlog::info("Performing micro-group pin cell calcuations");
  spdlog::info("Please wait...");
  set_logging_level(LogLevel::warn);

  // Preload all nuclides
  moderator_->load_nuclides(ndl_);
  for (const auto& pin : pins_)
    std::visit([this](const auto& P) { P->load_nuclides(ndl_); }, pin);

  // First, get all normal fuel pin indices (that don't need a buffer)
  std::vector<std::size_t> fp_inds;
  std::vector<std::size_t> other_inds;
  fp_inds.reserve(pins_.size());
  other_inds.reserve(pins_.size());
  auto needs_buffer = [](const auto& p) { return p->needs_buffer(); };
  for (std::size_t i = 0; i < pins_.size(); i++) {
    if (std::visit(needs_buffer, pins_[i])) {
      other_inds.push_back(i);
    } else {
      fp_inds.push_back(i);
    }
  }

  // Resize the cell and flux arrays
  pin_1d_cells.resize(pins_.size(), nullptr);
  pin_1d_fluxes.resize(pins_.size(), nullptr);

  // Now do all normal fuel pins
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

  // Now we compute all of the non-fuel cells, and fuel pins needing a buffer
  const double buffer_rad = std::sqrt(9. * pitch_ * pitch_ / PI);
#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(other_inds.size()); ii++) {
    const std::size_t i = other_inds[static_cast<std::size_t>(ii)];
    const double fuel_dancoff = fuel_dancoff_corrections_[i];
    const double clad_dancoff = clad_dancoff_corrections_[i];

    const auto& pin_var = pins_[i];
    if (std::holds_alternative<std::shared_ptr<FuelPin>>(pin_var)) {
      const auto& pin = std::get<std::shared_ptr<FuelPin>>(pin_var);
      pin_1d_cells[i] =
          pin->make_cylindrical_cell(pitch_, buffer_rad, avg_fp_, fuel_dancoff,
                                     moderator_xs_, ndl_, clad_dancoff);
    } else if (std::holds_alternative<std::shared_ptr<GuideTube>>(pin_var)) {
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

  pin_cell_calc_timer.stop();
  set_logging_level(LogLevel::info);
  spdlog::info("Pin cell calculation time: {:.5E} s",
               pin_cell_calc_timer.elapsed_time());
}

void PWRAssembly::condense_xs() {
  spdlog::info("");
  spdlog::info("Performing pin cell energy condensation");
  set_logging_level(LogLevel::warn);

  auto needs_buffer = [](const auto& p) { return p->needs_buffer(); };

  for (int ii = 0; ii < static_cast<int>(pins_.size()); ii++) {
    const std::size_t i = static_cast<std::size_t>(ii);
    auto& pin = pins_[i];
    const auto& cell_flux = pin_1d_fluxes[i];

    std::size_t NR = cell_flux->nregions();
    if (std::visit(needs_buffer, pin)) NR--;

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

void PWRAssembly::moc_calc() {
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
  moc_geom_ = std::make_shared<Cartesian2D>(dx, dy);
  moc_geom_->set_tiles(moc_pins);

  moc_ = std::make_shared<MOCDriver>(moc_geom_, boundary_conditions_,
                                     boundary_conditions_, boundary_conditions_,
                                     boundary_conditions_, anisotropic_);
  if (plot_assembly_) {
    ImApp::App guiplotter(1920, 1080, "Scarabee MOC Plotter");
    guiplotter.enable_docking();
    guiplotter.push_layer(std::make_unique<MOCPlotter>(moc_.get()));
    guiplotter.run();
  }
  moc_->generate_tracks(num_azimuthal_angles_, track_spacing_,
                        polar_quadrature_);
  moc_->set_keff_tolerance(keff_tolerance_);
  moc_->set_flux_tolerance(flux_tolerance_);
  moc_->solve();
}

void PWRAssembly::criticality_spectrum() {
  if (criticality_spectrum_method_.has_value() == false) {
    return;
  }

  spdlog::info("");
  spdlog::info("Performing {} criticality spectrum calculation",
               criticality_spectrum_method_.value());

  const auto homogenized_moc = moc_->homogenize();

  if (criticality_spectrum_method_.value() == "P1") {
    criticality_spectrum_ =
        std::make_shared<P1CriticalitySpectrum>(homogenized_moc);
  } else {
    criticality_spectrum_ =
        std::make_shared<B1CriticalitySpectrum>(homogenized_moc);
  }

  moc_->apply_criticality_spectrum(criticality_spectrum_->flux());

  spdlog::info("Kinf    : {:.5f}", criticality_spectrum_->k_inf());
  spdlog::info("Buckling: {:.5f}", criticality_spectrum_->buckling());
}

void PWRAssembly::compute_form_factors() {
  form_factors_ = xt::zeros<double>({shape_.first, shape_.second});

  const Direction u(1., 0.);

  for (std::size_t j = 0; j < shape_.second; j++) {
    const double y = moc_->y_max() - (static_cast<double>(j) + 0.5) * pitch_;
    for (std::size_t i = 0; i < shape_.first; i++) {
      const double x = moc_->x_min() + (static_cast<double>(i) + 0.5) * pitch_;
      const Vector r(x, y);
      auto inds = moc_->get_all_fsr_in_cell(r, u);

      // We now compute power in the cell
      for (auto fi : inds) {
        const double V = moc_->volume(fi);
        const auto& xs = moc_->xs(fi);
        for (std::size_t g = 0; g < xs->ngroups(); g++) {
          form_factors_(j, i) += V * xs->Ef(g) * moc_->flux(fi, g);
        }
      }
    }
  }

  // Now we compute the average pin power, and normalize the form factors
  const double avg_pwr = xt::mean(form_factors_)();
  form_factors_ /= avg_pwr;
}

void PWRAssembly::few_group_xs() {
  spdlog::info("");
  spdlog::info("Generating few group cross sections");

  const auto homog_xs = moc_->homogenize();
  const auto diff_xs = homog_xs->diffusion_xs();
  const auto flux_spectrum = moc_->homogenize_flux_spectrum();
  auto NG = homog_xs->ngroups();
  const bool fissile = homog_xs->fissile();

  xt::xtensor<double, 1> D = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 2> Es = xt::zeros<double>({NG, NG});
  xt::xtensor<double, 1> Ef, vEf, chi;

  // According to Smith, one should do energy condensation on the
  // diffusion coefficients, and not on the transport cross sections which
  // one could then use to make diffusion coefficients [1]. This is in
  // contradiction to Lattice Physics Computations which states that
  // either method is acceptable [2]. In light of these comments, I have
  // chosen to go with Smith's recommendation of performing energy
  // condensation on the diffusion coefficients.
  diffusion_xs_ =
      diff_xs->condense(few_group_condensation_scheme_, flux_spectrum);

  NG = diffusion_xs_->ngroups();

  D = xt::zeros<double>({NG});
  Ea = xt::zeros<double>({NG});
  Es = xt::zeros<double>({NG, NG});
  Ef = xt::zeros<double>({NG});
  vEf = xt::zeros<double>({NG});
  chi = xt::zeros<double>({NG});

  for (std::size_t g = 0; g < NG; g++) {
    D(g) = diffusion_xs_->D(g);
    Ea(g) = diffusion_xs_->Ea(g);
    Ef(g) = diffusion_xs_->Ef(g);
    vEf(g) = diffusion_xs_->vEf(g);
    chi(g) = diffusion_xs_->chi(g);

    for (std::size_t gg = 0; gg < NG; gg++) {
      Es(g, gg) = diffusion_xs_->Es(g, gg);
    }
  }

  std::stringstream D_str, Ea_str, Ef_str, vEf_str, chi_str, Es_str;
  D_str << "D  : " << D;
  Ea_str << "Ea : " << Ea;
  Ef_str << "Ef : " << Ef;
  vEf_str << "vEf: " << vEf;
  chi_str << "chi: " << chi;

  spdlog::info(D_str.str());
  spdlog::info(Ea_str.str());
  if (fissile) {
    spdlog::info(Ef_str.str());
    spdlog::info(vEf_str.str());
    spdlog::info(chi_str.str());
  }

  Es_str << "Es : " << xt::view(Es, 0, xt::all());
  spdlog::info(Es_str.str());
  for (std::size_t g = 1; g < NG; g++) {
    Es_str.str(std::string());
    Es_str << "     " << xt::view(Es, g, xt::all());
    spdlog::info(Es_str.str());
  }
}

void PWRAssembly::compute_adf_cdf() {
  // Number of few groups
  const std::size_t NG = few_group_condensation_scheme_.size();

  // We first need to compute the homogenous few-group flux for the assembly
  std::vector<double> homog_flx(NG, 0.);
  double total_volume = 0.;
  for (std::size_t i = 0; i < moc_->nfsr(); i++) {
    const double Vi = moc_->volume(i);
    total_volume += Vi;

    for (std::size_t G = 0; G < NG; G++) {
      const std::size_t gmin = few_group_condensation_scheme_[G].first;
      const std::size_t gmax = few_group_condensation_scheme_[G].first;

      for (std::size_t g = gmin; g <= gmax; g++) {
        homog_flx[G] += Vi * moc_->flux(i, g);
      }
    }
  }
  for (auto& flx : homog_flx) flx /= total_volume;

  // Trace the needed segments along the assembly
  const auto xp_segments = moc_->trace_fsr_segments(
      Vector(moc_->x_max() - 0.01, moc_->y_max()), Direction(0., -1.));

  const auto xn_segments = moc_->trace_fsr_segments(
      Vector(moc_->x_min() + 0.01, moc_->y_max()), Direction(0., -1.));

  const auto yp_segments = moc_->trace_fsr_segments(
      Vector(moc_->x_min(), moc_->y_max() - 0.01), Direction(1., 0.));

  const auto yn_segments = moc_->trace_fsr_segments(
      Vector(moc_->x_min(), moc_->y_min() + 0.01), Direction(1., 0.));

  // Get average flux in each macro group along each surface
  const auto xp_flx = compute_avg_surface_flx(xp_segments);
  const auto xn_flx = compute_avg_surface_flx(xn_segments);
  const auto yp_flx = compute_avg_surface_flx(yp_segments);
  const auto yn_flx = compute_avg_surface_flx(yn_segments);

  // Load the adf array
  adf_ = xt::zeros<double>({NG, static_cast<std::size_t>(4)});
  for (std::size_t G = 0; G < NG; G++) {
    adf_(G, DiffusionData::ADF::YP) = yp_flx[G] / homog_flx[G];
    adf_(G, DiffusionData::ADF::XP) = xp_flx[G] / homog_flx[G];
    adf_(G, DiffusionData::ADF::YN) = yn_flx[G] / homog_flx[G];
    adf_(G, DiffusionData::ADF::XN) = xn_flx[G] / homog_flx[G];
  }

  // We now need to compute the corner fluxes for the cdf
  const auto I_flx = compute_avg_flx(
      Vector(moc_->x_max() - 0.01, moc_->y_max() - 0.011), Direction(-1., -1.));
  const auto II_flx = compute_avg_flx(
      Vector(moc_->x_min() + 0.01, moc_->y_max() - 0.011), Direction(1., -1.));
  const auto III_flx = compute_avg_flx(
      Vector(moc_->x_min() + 0.01, moc_->y_min() + 0.011), Direction(1., 1.));
  const auto IV_flx = compute_avg_flx(
      Vector(moc_->x_max() - 0.01, moc_->y_min() + 0.011), Direction(-1., 1.));
  cdf_ = xt::zeros<double>({NG, static_cast<std::size_t>(4)});
  for (std::size_t G = 0; G < NG; G++) {
    cdf_(G, DiffusionData::CDF::I) = I_flx[G] / homog_flx[G];
    cdf_(G, DiffusionData::CDF::II) = II_flx[G] / homog_flx[G];
    cdf_(G, DiffusionData::CDF::III) = III_flx[G] / homog_flx[G];
    cdf_(G, DiffusionData::CDF::IV) = IV_flx[G] / homog_flx[G];
  }

  std::stringstream str;
  str << "ADF : " << xt::view(adf_, 0, xt::all());
  spdlog::info(str.str());
  for (std::size_t G = 1; G < NG; G++) {
    str.str(std::string());
    str << "      " << xt::view(adf_, G, xt::all());
    spdlog::info(str.str());
  }
  str.str(std::string());

  str << "CDF : " << xt::view(cdf_, 0, xt::all());
  spdlog::info(str.str());
  for (std::size_t G = 1; G < NG; G++) {
    str.str(std::string());
    str << "      " << xt::view(cdf_, G, xt::all());
    spdlog::info(str.str());
  }
}

void PWRAssembly::save_diffusion_data(const std::string& fname) const {
  if (diffusion_xs_ == nullptr || form_factors_.size() == 0 ||
      adf_.size() == 0 || cdf_.size() == 0) {
    auto mssg = "Cannot save DiffusionData. Assembly has not been solved.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  DiffusionData dd(diffusion_xs_, form_factors_, adf_, cdf_);
  dd.save(fname);
}

std::vector<double> PWRAssembly::compute_avg_surface_flx(
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

std::vector<double> PWRAssembly::compute_avg_flx(const Vector& r,
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

void PWRAssembly::save(const std::string& fname) const {
  if (std::filesystem::exists(fname)) {
    std::filesystem::remove(fname);
  }

  std::ofstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryOutputArchive arc(file);

  arc(*this);
}

std::shared_ptr<PWRAssembly> PWRAssembly::load(const std::string& fname) {
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  std::shared_ptr<PWRAssembly> out(new PWRAssembly());

  std::ifstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryInputArchive arc(file);

  arc(*out);

  return out;
}

}  // namespace scarabee

// REFERENCES
// [1] K. S. Smith, “Nodal diffusion methods and lattice physics data in LWR
//     analyses: Understanding numerous subtle details,” Prog Nucl Energ,
//     vol. 101, pp. 360–369, 2017, doi: 10.1016/j.pnucene.2017.06.013.
// [2] D. Knott and A. Yamamoto, "Lattice Physics Computations" in
//     Handbook of Nuclear Engineering, 2010, p 1226.
