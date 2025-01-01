#include <data/nd_library.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <cstdlib>
#include <optional>
#include <sstream>

namespace scarabee {

void NuclideHandle::load_xs_from_hdf5(const NDLibrary& ndl, std::size_t max_l) {
  if (this->loaded()) return;

  auto grp = ndl.h5()->getGroup(this->name);

  // Create and allocate arrays
  absorption = std::make_shared<xt::xtensor<double, 3>>();
  absorption->resize({temperatures.size(), dilutions.size(), ndl.ngroups()});

  transport_correction = std::make_shared<xt::xtensor<double, 3>>();
  transport_correction->resize(
      {temperatures.size(), dilutions.size(), ndl.ngroups()});

  // Must read in packing now to know full length of arrays
  packing = std::make_shared<xt::xtensor<std::uint32_t, 2>>();
  packing->resize({ndl.ngroups(), static_cast<std::size_t>(3)});
  grp.getDataSet("matrix-compression").read_raw<std::uint32_t>(packing->data());
  const std::size_t len_scat_data = (*packing)(ndl.ngroups() - 1, 0) +
                                    (*packing)(ndl.ngroups() - 1, 2) + 1 -
                                    (*packing)(ndl.ngroups() - 1, 1);

  scatter = std::make_shared<xt::xtensor<double, 3>>();
  scatter->resize({temperatures.size(), dilutions.size(), len_scat_data});

  if (grp.exist("p1-scatter") && max_l >= 1) {
    p1_scatter = std::make_shared<xt::xtensor<double, 3>>();
    p1_scatter->resize({temperatures.size(), dilutions.size(), len_scat_data});
  }
  if (grp.exist("p2-scatter") && max_l >= 2) {
    p2_scatter = std::make_shared<xt::xtensor<double, 3>>();
    p2_scatter->resize({temperatures.size(), dilutions.size(), len_scat_data});
  }
  if (grp.exist("p3-scatter") && max_l >= 3) {
    p3_scatter = std::make_shared<xt::xtensor<double, 3>>();
    p3_scatter->resize({temperatures.size(), dilutions.size(), len_scat_data});
  }

  if (this->fissile) {
    fission = std::make_shared<xt::xtensor<double, 3>>();
    fission->resize({temperatures.size(), dilutions.size(), ndl.ngroups()});
    nu = std::make_shared<xt::xtensor<double, 1>>();
    nu->resize({ndl.ngroups()});
    chi = std::make_shared<xt::xtensor<double, 1>>();
    chi->resize({ndl.ngroups()});
  }

  // Read in data
  grp.getDataSet("absorption").read_raw<double>(absorption->data());
  grp.getDataSet("transport-correction")
      .read_raw<double>(transport_correction->data());
  grp.getDataSet("scatter").read_raw<double>(scatter->data());
  if (grp.exist("p1-scatter") && max_l >= 1) {
    grp.getDataSet("p1-scatter").read_raw<double>(p1_scatter->data());
  }
  if (grp.exist("p2-scatter") && max_l >= 2) {
    grp.getDataSet("p2-scatter").read_raw<double>(p2_scatter->data());
  }
  if (grp.exist("p3-scatter") && max_l >= 3) {
    grp.getDataSet("p3-scatter").read_raw<double>(p3_scatter->data());
  }
  if (this->fissile) {
    grp.getDataSet("fission").read_raw<double>(fission->data());
    grp.getDataSet("nu").read_raw<double>(nu->data());
    grp.getDataSet("chi").read_raw<double>(chi->data());
  }
}

void NuclideHandle::unload() {
  absorption = nullptr;
  transport_correction = nullptr;
  scatter = nullptr;
  p1_scatter = nullptr;
  p2_scatter = nullptr;
  p3_scatter = nullptr;
  fission = nullptr;
  chi = nullptr;
  nu = nullptr;
}

NDLibrary::NDLibrary()
    : nuclide_handles_(),
      group_bounds_(),
      macro_group_condensation_scheme_(std::nullopt),
      few_group_condensation_scheme_(std::nullopt),
      library_(),
      group_structure_(),
      ngroups_(0),
      h5_(nullptr) {
  // Get the environment variable
  const char* ndl_env = std::getenv(NDL_ENV_VAR);
  if (ndl_env == nullptr) {
    auto mssg = "Environment variable " NDL_ENV_VAR " is not set.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Get the string for the file name.
  std::string fname(ndl_env);

  // Make sure HDF5 file exists
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Open the HDF5 file
  h5_ = std::make_shared<H5::File>(fname, H5::File::ReadOnly);

  this->init();
}

NDLibrary::NDLibrary(const std::string& fname)
    : nuclide_handles_(),
      group_bounds_(),
      macro_group_condensation_scheme_(std::nullopt),
      few_group_condensation_scheme_(std::nullopt),
      library_(),
      group_structure_(),
      ngroups_(0),
      h5_(nullptr) {
  // Make sure HDF5 file exists
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Open the HDF5 file
  h5_ = std::make_shared<H5::File>(fname, H5::File::ReadOnly);

  this->init();
}

void NDLibrary::init() {
  // Get info on library
  if (h5_->hasAttribute("library"))
    library_ = h5_->getAttribute("library").read<std::string>();
  if (h5_->hasAttribute("group-structure"))
    group_structure_ = h5_->getAttribute("group-structure").read<std::string>();
  if (h5_->hasAttribute("group-bounds"))
    group_bounds_ =
        h5_->getAttribute("group-bounds").read<std::vector<double>>();
  if (h5_->hasAttribute("ngroups"))
    ngroups_ = h5_->getAttribute("ngroups").read<std::size_t>();

  // Check if default condensation schemes are provided
  auto get_cond_scheme =
      [this](const std::string& key,
             std::optional<std::vector<std::pair<std::size_t, std::size_t>>>&
                 scheme,
             const std::string& mssg) {
        if (h5_->hasAttribute(key)) {
          const auto attr = h5_->getAttribute(key);
          const auto dims = attr.getMemSpace().getDimensions();

          if (dims.size() != 2 || dims[1] != 2) {
            spdlog::error(mssg);
            throw ScarabeeException(mssg);
          }

          xt::xtensor<double, 2> cond_scheme =
              xt::zeros<double>({dims[0], dims[1]});
          attr.read_raw<double>(cond_scheme.data());

          scheme =
              std::vector<std::pair<std::size_t, std::size_t>>(dims[0], {0, 0});
          for (std::size_t G = 0; G < dims[0]; G++) {
            (*scheme)[G].first = static_cast<std::size_t>(cond_scheme(G, 0));
            (*scheme)[G].second = static_cast<std::size_t>(cond_scheme(G, 1));
          }
        }
      };

  get_cond_scheme("macro-group-condensation-scheme",
                  macro_group_condensation_scheme_,
                  "Nuclear data library provided macro-group condensation "
                  "scheme has an invalide shape.");
  get_cond_scheme("few-group-condensation-scheme",
                  few_group_condensation_scheme_,
                  "Nuclear data library provided few-group condensation scheme "
                  "has an invalide shape.");
  get_cond_scheme("reflector-few-group-condensation-scheme",
                  reflector_few_group_condensation_scheme_,
                  "Nuclear data library provided reflector few-group "
                  "condensation scheme has an invalide shape.");

  // Read all nuclide handles
  auto nuc_names = h5_->listObjectNames();
  for (const auto& nuc : nuc_names) {
    auto grp = h5_->getGroup(nuc);

    nuclide_handles_.emplace(std::make_pair(nuc, NuclideHandle()));
    auto& handle = nuclide_handles_.at(nuc);
    handle.name = nuc;

    // Read nuclide info
    handle.label = grp.getAttribute("label").read<std::string>();
    handle.temperatures =
        grp.getAttribute("temperatures").read<std::vector<double>>();
    handle.awr = grp.getAttribute("awr").read<double>();
    handle.potential_xs = grp.getAttribute("potential-xs").read<double>();

    // Intermediate resonance parameter
    if (grp.hasAttribute("ir-lambda")) {
      handle.ir_lambda = grp.getAttribute("ir-lambda").read<double>();
    } else {
      handle.ir_lambda = 1.;
    }

    handle.ZA = grp.getAttribute("ZA").read<std::uint32_t>();
    handle.fissile = grp.getAttribute("fissile").read<bool>();
    handle.resonant = grp.getAttribute("resonant").read<bool>();
    handle.dilutions =
        grp.getAttribute("dilutions").read<std::vector<double>>();
  }
}

const NuclideHandle& NDLibrary::get_nuclide(const std::string& name) const {
  if (nuclide_handles_.find(name) == nuclide_handles_.end()) {
    std::stringstream mssg;
    mssg << "Could not find nuclde by name of \"" << name << "\".";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return nuclide_handles_.at(name);
}

void NDLibrary::unload() {
  for (auto& nuc_handle : nuclide_handles_) {
    nuc_handle.second.unload();
  }
}

NuclideHandle& NDLibrary::get_nuclide(const std::string& name) {
  if (nuclide_handles_.find(name) == nuclide_handles_.end()) {
    std::stringstream mssg;
    mssg << "Could not find nuclde by name of \"" << name << "\".";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return nuclide_handles_.at(name);
}

std::shared_ptr<CrossSection> NDLibrary::interp_xs(const std::string& name,
                                                   const double temp,
                                                   const double dil,
                                                   std::size_t max_l) {
  auto& nuc = this->get_nuclide(name);

  // Get temperature interpolation factors
  std::size_t it = 0;  // temperature index
  double f_temp = 0.;  // temperature interpolation factor
  get_temp_interp_params(temp, nuc, it, f_temp);

  // Get dilution interpolation factors
  std::size_t id = 0;  // dilution index
  double f_dil = 0.;   // dilution interpolation factor
  get_dil_interp_params(dil, nuc, id, f_dil);

  if (nuc.loaded() == false) {
    nuc.load_xs_from_hdf5(*this, max_l);
  }

  //--------------------------------------------------------
  // Do transport correction
  xt::xtensor<double, 1> Dtr;
  this->interp_1d(Dtr, *nuc.transport_correction, it, f_temp, id, f_dil);

  //--------------------------------------------------------
  // Do absorption interpolation
  xt::xtensor<double, 1> Ea;
  this->interp_1d(Ea, *nuc.absorption, it, f_temp, id, f_dil);

  //--------------------------------------------------------
  // Create the scattering matrices
  if (max_l == 3 && nuc.p3_scatter == nullptr) max_l--;
  if (max_l == 2 && nuc.p2_scatter == nullptr) max_l--;
  if (max_l == 1 && nuc.p1_scatter == nullptr) max_l--;
  xt::xtensor<double, 2> Es =
      xt::zeros<double>({max_l + 1, nuc.scatter->shape()[2]});

  //--------------------------------------------------------
  // Do P0 scattering interpolation
  xt::xtensor<double, 1> temp_EsPl;
  this->interp_1d(temp_EsPl, *nuc.scatter, it, f_temp, id, f_dil);
  xt::view(Es, 0, xt::all()) = temp_EsPl;

  //--------------------------------------------------------
  // Do P1 scattering interpolation
  if (nuc.p1_scatter) {
    this->interp_1d(temp_EsPl, *nuc.p1_scatter, it, f_temp, id, f_dil);
    xt::view(Es, 1, xt::all()) = temp_EsPl;
  }

  //--------------------------------------------------------
  // Do P2 scattering interpolation
  if (nuc.p2_scatter && max_l >= 2) {
    this->interp_1d(temp_EsPl, *nuc.p2_scatter, it, f_temp, id, f_dil);
    xt::view(Es, 2, xt::all()) = temp_EsPl;
  }

  //--------------------------------------------------------
  // Do P3 scattering interpolation
  if (nuc.p3_scatter && max_l >= 3) {
    this->interp_1d(temp_EsPl, *nuc.p3_scatter, it, f_temp, id, f_dil);
    xt::view(Es, 3, xt::all()) = temp_EsPl;
  }
  XS2D Es_xs2d(Es, *nuc.packing);

  //--------------------------------------------------------
  // Do fission interpolation
  xt::xtensor<double, 1> Ef = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> nu = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> chi = xt::zeros<double>({ngroups_});
  if (nuc.fissile) {
    this->interp_1d(Ef, *nuc.fission, it, f_temp, id, f_dil);
    nu = *nuc.nu;
    chi = *nuc.chi;
  }

  // Reconstruct total
  xt::xtensor<double, 1> Et = xt::zeros<double>({ngroups_});
  for (std::size_t g = 0; g < ngroups_; g++) {
    Et(g) = Ea(g) + Es_xs2d(0, g);
  }

  // Make temp CrossSection
  return std::make_shared<CrossSection>(XS1D(Et), XS1D(Dtr), XS1D(Ea), Es_xs2d,
                                        XS1D(Ef), XS1D(nu * Ef), XS1D(chi));
}

std::shared_ptr<CrossSection> NDLibrary::two_term_xs(
    const std::string& name, const double temp, const double b1,
    const double b2, const double bg_xs_1, const double bg_xs_2,
    std::size_t max_l) {
  // See reference [1] to understand this interpolation scheme, in addition to
  // the calculation of the flux based on the pot_xs and sig_a.

  // Get the two cross section sets
  auto xs_1 = interp_xs(name, temp, bg_xs_1, max_l);
  auto xs_2 = interp_xs(name, temp, bg_xs_2, max_l);

  const auto& nuclide = get_nuclide(name);
  const double pot_xs = nuclide.ir_lambda * nuclide.potential_xs;

  XS1D Et(xt::zeros<double>({ngroups_}));
  XS1D Dtr(xt::zeros<double>({ngroups_}));
  XS1D Ea(xt::zeros<double>({ngroups_}));
  XS2D Es = xs_1->Es_XS2D().zeros_like();
  XS1D Ef(xt::zeros<double>({ngroups_}));
  XS1D vEf(xt::zeros<double>({ngroups_}));
  XS1D chi(xt::zeros<double>({ngroups_}));

  for (std::size_t g = 0; g < ngroups_; g++) {
    // Calculate the two flux values
    const double flux_1_g =
        (pot_xs + bg_xs_1) / (xs_1->Ea(g) + pot_xs + bg_xs_1);
    const double flux_2_g =
        (pot_xs + bg_xs_2) / (xs_2->Ea(g) + pot_xs + bg_xs_2);

    // Calcualte the two weighting factors
    const double f1_g = b1 * flux_1_g / (b1 * flux_1_g + b2 * flux_2_g);
    const double f2_g = b2 * flux_2_g / (b1 * flux_1_g + b2 * flux_2_g);

    // Compute the xs values
    Dtr.set_value(g, f1_g * xs_1->Dtr(g) + f2_g * xs_2->Dtr(g));
    Ea.set_value(g, f1_g * xs_1->Ea(g) + f2_g * xs_2->Ea(g));
    Ef.set_value(g, f1_g * xs_1->Ef(g) + f2_g * xs_2->Ef(g));

    // Min and Max outgoing groups
    const std::size_t gg_min = static_cast<std::size_t>(Es.packing()(g, 1));
    const std::size_t gg_max = static_cast<std::size_t>(Es.packing()(g, 2));
    for (std::size_t l = 0; l <= xs_1->max_legendre_order(); l++) {
      for (std::size_t g_out = gg_min; g_out <= gg_max; g_out++) {
        Es.set_value(
            l, g, g_out,
            f1_g * xs_1->Es(l, g, g_out) + f2_g * xs_2->Es(l, g, g_out));
      }
    }
    Et.set_value(g, Ea(g) + Es(0, g));

    const double vEf1 = f1_g * xs_1->vEf(g);
    const double vEf2 = f2_g * xs_2->vEf(g);
    vEf.set_value(g, vEf1 + vEf2);

    // Chi isn't stored on temp or dilution, so just assign value
    chi.set_value(g, xs_1->chi(g));
  }

  return std::make_shared<CrossSection>(Et, Dtr, Ea, Es, Ef, vEf, chi);
}

std::shared_ptr<CrossSection> NDLibrary::ring_two_term_xs(
    const std::string& name, const double temp, const double a1,
    const double a2, const double b1, const double b2, const double mat_pot_xs,
    const double N, const double Rfuel, const double Rin, const double Rout,
    std::size_t max_l) {
  if (Rin >= Rout) {
    auto mssg = "Rin must be < Rout.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Rout > Rfuel) {
    auto mssg = "Rout must be < Rfuel.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const auto& nuclide = get_nuclide(name);
  const double pot_xs = nuclide.ir_lambda * nuclide.potential_xs;
  const double macro_pot_xs = N * pot_xs;

  if (max_l == 3 && nuclide.p3_scatter == nullptr) max_l--;
  if (max_l == 2 && nuclide.p2_scatter == nullptr) max_l--;
  if (max_l == 1 && nuclide.p1_scatter == nullptr) max_l--;

  XS1D Et(xt::zeros<double>({ngroups_}));
  XS1D Dtr(xt::zeros<double>({ngroups_}));
  XS1D Ea(xt::zeros<double>({ngroups_}));
  std::optional<XS2D> Es{std::nullopt};
  XS1D Ef(xt::zeros<double>({ngroups_}));
  XS1D vEf(xt::zeros<double>({ngroups_}));
  XS1D chi(xt::zeros<double>({ngroups_}));
  if (nuclide.chi) {
    for (std::size_t g = 0; g < ngroups_; g++) {
      chi.set_value(g, (*nuclide.chi)(g));
    }
  }

  // Denominators of the weighting factor for each energy group.
  xt::xtensor<double, 1> denoms = xt::zeros<double>({ngroups_});

  for (std::size_t m = 1; m <= 4; m++) {
    const std::pair<double, double> eta_lm = this->eta_lm(m, Rfuel, Rin, Rout);
    const double eta_m = eta_lm.first;
    const double l_m = eta_lm.second;

    // If eta_m is zero, then this m has no contribution to the xs.
    // This happens for the last ring of a pin.
    if (eta_m == 0.) continue;

    // Calculate the background xs
    const double bg_xs_1 =
        l_m > 0. ? (mat_pot_xs - macro_pot_xs + a1 / l_m) / N : 1.E10;
    const double bg_xs_2 =
        l_m > 0. ? (mat_pot_xs - macro_pot_xs + a2 / l_m) / N : 1.E10;

    // Get the two cross section sets
    auto xs_1 = interp_xs(name, temp, bg_xs_1, max_l);
    auto xs_2 = interp_xs(name, temp, bg_xs_2, max_l);

    // Now that we have a xs instance, initialize the scatter matrix
    if (Es.has_value() == false) Es = xs_1->Es_XS2D().zeros_like();

    for (std::size_t g = 0; g < ngroups_; g++) {
      // Calculate the two flux values
      const double flux_1_g =
          (pot_xs + bg_xs_1) / (xs_1->Ea(g) + pot_xs + bg_xs_1);
      const double flux_2_g =
          (pot_xs + bg_xs_2) / (xs_2->Ea(g) + pot_xs + bg_xs_2);

      // Add contributions to the denominator
      denoms(g) += eta_m * (b1 * flux_1_g + b2 * flux_2_g);

      // Add contributions to the xs
      // Compute the xs values
      Dtr.set_value(g, Dtr(g) + eta_m * (b1 * flux_1_g * xs_1->Dtr(g) +
                                         b2 * flux_2_g * xs_2->Dtr(g)));
      Ea.set_value(g, Ea(g) + eta_m * (b1 * flux_1_g * xs_1->Ea(g) +
                                       b2 * flux_2_g * xs_2->Ea(g)));
      Ef.set_value(g, Ef(g) + eta_m * (b1 * flux_1_g * xs_1->Ef(g) +
                                       b2 * flux_2_g * xs_2->Ef(g)));
      vEf.set_value(g, vEf(g) + eta_m * (b1 * flux_1_g * xs_1->vEf(g) +
                                         b2 * flux_2_g * xs_2->vEf(g)));
      for (std::size_t l = 0; l <= max_l; l++) {
        for (std::size_t g_out = 0; g_out < ngroups_; g_out++) {
          const double new_val =
              (*Es)(l, g, g_out) +
              eta_m * (b1 * flux_1_g * xs_1->Es(l, g, g_out) +
                       b2 * flux_2_g * xs_2->Es(l, g, g_out));
          if (new_val != 0.) Es->set_value(l, g, g_out, new_val);
        }  // For all outgoing groups
      }
    }  // For all groups
  }    // For 4 lumps

  // Now we go through and normalize each group by the denom, and calculate Et
  for (std::size_t g = 0; g < ngroups_; g++) {
    const double invs_denom = 1. / denoms(g);
    Dtr.set_value(g, Dtr(g) * invs_denom);
    Ea.set_value(g, Ea(g) * invs_denom);
    Ef.set_value(g, Ef(g) * invs_denom);
    vEf.set_value(g, vEf(g) * invs_denom);

    for (std::size_t l = 0; l <= max_l; l++) {
      for (std::size_t g_out = 0; g_out < ngroups_; g_out++) {
        auto Es_g_gout = (*Es)(l, g, g_out);
        if (Es_g_gout != 0.) Es->set_value(l, g, g_out, Es_g_gout * invs_denom);
      }
    }

    Et.set_value(g, Ea(g) + (*Es)(0, g));
  }

  return std::make_shared<CrossSection>(Et, Dtr, Ea, *Es, Ef, vEf, chi);
}

void NDLibrary::get_temp_interp_params(double temp, const NuclideHandle& nuc,
                                       std::size_t& i, double& f) const {
  if (temp <= nuc.temperatures.front()) {
    i = 0;
    f = 0.;
    return;
  } else if (temp >= nuc.temperatures.back()) {
    i = nuc.temperatures.size() - 2;
    f = 1.;
    return;
  }

  for (i = 0; i < nuc.temperatures.size() - 1; i++) {
    double T_i = nuc.temperatures[i];
    double T_i1 = nuc.temperatures[i + 1];

    if (temp >= T_i && temp <= T_i1) {
      f = (std::sqrt(temp) - std::sqrt(T_i)) /
          (std::sqrt(T_i1) - std::sqrt(T_i));
      break;
    }
  }
  if (f < 0.)
    f = 0.;
  else if (f > 1.)
    f = 1.;
}

void NDLibrary::get_dil_interp_params(double dil, const NuclideHandle& nuc,
                                      std::size_t& i, double& f) const {
  if (dil <= nuc.dilutions.front()) {
    i = 0;
    f = 0.;
    return;
  } else if (dil >= nuc.dilutions.back()) {
    i = nuc.dilutions.size() - 2;
    f = 1.;
    return;
  }

  for (i = 0; i < nuc.dilutions.size() - 1; i++) {
    double d_i = nuc.dilutions[i];
    double d_i1 = nuc.dilutions[i + 1];

    if (dil >= d_i && dil <= d_i1) {
      f = (dil - d_i) / (d_i1 - d_i);
      break;
    }
  }
  if (f < 0.)
    f = 0.;
  else if (f > 1.)
    f = 1.;
}

void NDLibrary::interp_1d(xt::xtensor<double, 1>& E,
                          const xt::xtensor<double, 2>& nE, std::size_t it,
                          double f_temp) const {
  if (f_temp > 0.) {
    E = (1. - f_temp) * xt::view(nE, it, xt::all()) +
        f_temp * xt::view(nE, it + 1, xt::all());
  } else {
    E = xt::view(nE, it, xt::all());
  }
}

void NDLibrary::interp_1d(xt::xtensor<double, 1>& E,
                          const xt::xtensor<double, 3>& nE, std::size_t it,
                          double f_temp, std::size_t id, double f_dil) const {
  if (f_temp > 0.) {
    if (f_dil > 0.) {
      E = (1. - f_temp) * ((1. - f_dil) * xt::view(nE, it, id, xt::all()) +
                           f_dil * xt::view(nE, it, id + 1, xt::all())) +
          f_temp * ((1. - f_dil) * xt::view(nE, it + 1, id, xt::all()) +
                    f_dil * xt::view(nE, it + 1, id + 1, xt::all()));
    } else {
      E = (1. - f_temp) * xt::view(nE, it, id, xt::all()) +
          f_temp * xt::view(nE, it + 1, id, xt::all());
    }
  } else {
    if (f_dil > 0.) {
      E = (1. - f_dil) * xt::view(nE, it, id, xt::all()) +
          f_dil * xt::view(nE, it, id + 1, xt::all());
    } else {
      E = xt::view(nE, it, id, xt::all());
    }
  }
}

void NDLibrary::interp_2d(xt::xtensor<double, 2>& E,
                          const xt::xtensor<double, 4>& nE, std::size_t it,
                          double f_temp, std::size_t id, double f_dil) const {
  if (f_temp > 0.) {
    if (f_dil > 0.) {
      E = (1. - f_temp) *
              ((1. - f_dil) * xt::view(nE, it, id, xt::all(), xt::all()) +
               f_dil * xt::view(nE, it, id + 1, xt::all(), xt::all())) +
          f_temp *
              ((1. - f_dil) * xt::view(nE, it + 1, id, xt::all(), xt::all()) +
               f_dil * xt::view(nE, it + 1, id + 1, xt::all(), xt::all()));
    } else {
      E = (1. - f_temp) * xt::view(nE, it, id, xt::all(), xt::all()) +
          f_temp * xt::view(nE, it + 1, id, xt::all(), xt::all());
    }
  } else {
    if (f_dil > 0.) {
      E = (1. - f_dil) * xt::view(nE, it, id, xt::all(), xt::all()) +
          f_dil * xt::view(nE, it, id + 1, xt::all(), xt::all());
    } else {
      E = xt::view(nE, it, id, xt::all(), xt::all());
    }
  }
}

std::pair<double, double> NDLibrary::eta_lm(std::size_t m, double Rfuel,
                                            double Rin, double Rout) const {
  if (m == 0 || m > 4) {
    auto mssg = "Invalid m.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (Rin >= Rout) {
    auto mssg = "Rin >= Rout.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Shouldn't need to check m, as this is a private method
  const double p_i = std::min(Rout / Rfuel, 1.);
  const double p_im = Rin / Rfuel;

  double p = p_i;
  if (m == 3 || m == 4) p = p_im;

  double theta = 0.5 * PI * p;
  if (m == 2 || m == 4) theta = -theta;

  // l = 4V_ring / S_pin = 4 pi (Rout^2 - Rin^2) / (2 pi Rfuel)
  const double l = 2. * (Rout * Rout - Rin * Rin) / Rfuel;

  const double T1 = std::sqrt(1. - p * p);
  const double T2 = Rin > 0. ? std::asin(p) / p : 1.;

  const double lm = (2. * Rfuel / PI) * (T1 + T2 + theta);

  double eta = p * lm / l;
  if (m == 2 || m == 3) eta = -eta;

  return {eta, lm};
}

}  // namespace scarabee

// References
// [1] H. Koike, K. Yamaji, K. Kirimura, D. Sato, H. Matsumoto, and A. Yamamoto,
//     “Advanced resonance self-shielding method for gray resonance treatment in
//     lattice physics code GALAXY,” J. Nucl. Sci. Technol., vol. 49, no. 7,
//     pp. 725–747, 2012, doi: 10.1080/00223131.2012.693885.
