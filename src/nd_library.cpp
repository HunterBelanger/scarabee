#include <data/nd_library.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <filesystem>
#include <sstream>

namespace scarabee {

void NuclideHandle::load_xs_from_hdf5(const NDLibrary& ndl) {
  if (this->loaded()) return;

  auto grp = ndl.h5()->getGroup(this->name);

  // Create and allocate arrays
  absorption = std::make_shared<xt::xtensor<double, 3>>();
  absorption->resize({temperatures.size(), dilutions.size(), ndl.ngroups()});
  scatter = std::make_shared<xt::xtensor<double, 4>>();
  scatter->resize(
      {temperatures.size(), dilutions.size(), ndl.ngroups(), ndl.ngroups()});
  p1_scatter = std::make_shared<xt::xtensor<double, 4>>();
  p1_scatter->resize(
      {temperatures.size(), dilutions.size(), ndl.ngroups(), ndl.ngroups()});
  if (this->fissile) {
    fission = std::make_shared<xt::xtensor<double, 3>>();
    fission->resize({temperatures.size(), dilutions.size(), ndl.ngroups()});
    nu = std::make_shared<xt::xtensor<double, 2>>();
    nu->resize({temperatures.size(), ndl.ngroups()});
    chi = std::make_shared<xt::xtensor<double, 2>>();
    chi->resize({temperatures.size(), ndl.ngroups()});
  }

  // Read in data
  grp.getDataSet("absorption").read<double>(absorption->data());
  grp.getDataSet("scatter").read<double>(scatter->data());
  grp.getDataSet("p1-scatter").read<double>(p1_scatter->data());
  if (this->fissile) {
    grp.getDataSet("fission").read<double>(fission->data());
    grp.getDataSet("nu").read<double>(nu->data());
    grp.getDataSet("chi").read<double>(chi->data());
  }
}

NDLibrary::NDLibrary(const std::string& fname)
    : nuclide_handles_(),
      group_bounds_(),
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

NuclideHandle& NDLibrary::get_nuclide(const std::string& name) {
  if (nuclide_handles_.find(name) == nuclide_handles_.end()) {
    std::stringstream mssg;
    mssg << "Could not find nuclde by name of \"" << name << "\".";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return nuclide_handles_.at(name);
}

std::shared_ptr<CrossSection> NDLibrary::interp_nuclide_xs(
    const std::string& name, const double temp, const double dil) {
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
    nuc.load_xs_from_hdf5(*this);
  }

  //--------------------------------------------------------
  // Do absorption interpolation
  xt::xtensor<double, 1> Ea;
  this->interp_1d(Ea, *nuc.absorption, it, f_temp, id, f_dil);

  //--------------------------------------------------------
  // Do scattering interpolation
  xt::xtensor<double, 2> Es;
  this->interp_2d(Es, *nuc.scatter, it, f_temp, id, f_dil);

  //--------------------------------------------------------
  // Do p1 scattering interpolation
  xt::xtensor<double, 2> Es1;
  this->interp_2d(Es1, *nuc.p1_scatter, it, f_temp, id, f_dil);

  //--------------------------------------------------------
  // Do fission interpolation
  xt::xtensor<double, 1> Ef = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> nu = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> chi = xt::zeros<double>({ngroups_});
  if (nuc.fissile) {
    this->interp_1d(Ef, *nuc.fission, it, f_temp, id, f_dil);
    this->interp_1d(nu, *nuc.nu, it, f_temp);
    this->interp_1d(chi, *nuc.chi, it, f_temp);
  }

  // Reconstruct total, removing p1
  xt::xtensor<double, 1> Et = xt::zeros<double>({ngroups_});
  for (std::size_t g = 0; g < ngroups_; g++) {
    Et(g) = Ea(g) + xt::sum(xt::view(Es, g, xt::all()))();
    Et(g) -= Es1(g, g);
    Es(g, g) -= Es1(g, g);
  }

  // Make temp CrossSection
  return std::make_shared<CrossSection>(Et, Ea, Es, Ef, nu * Ef, chi);
}

std::shared_ptr<CrossSection> NDLibrary::carlvik_two_term(
    const std::string& name, const double mat_pot_xs, const double temp,
    const double N, const double C, const double Ee) {
  // This implementation is based on the methods outlined by Koike and Gibson
  // in his PhD thesis [1,2]. We start by computing the coefficients for the
  // two-term rational approximation, modified according to the Dancoff
  // correction factor, C.
  const double a1 = 0.5 * (C + 5. - std::sqrt(C * C + 34. * C + 1.));
  const double a2 = 0.5 * (C + 5. + std::sqrt(C * C + 34. * C + 1.));
  const double b1 = (a2 - (1. - C)) / (a2 - a1);
  const double b2 = 1. - b1;

  const double pot_xs = get_nuclide(name).potential_xs;
  const double macro_pot_xs = N * pot_xs;

  // Compute the two background xs values for the nuclide
  const double bg_xs_1 = (mat_pot_xs - macro_pot_xs + a1 * Ee) / N;
  const double bg_xs_2 = (mat_pot_xs - macro_pot_xs + a2 * Ee) / N;

  // Get the two cross section sets
  auto xs_1 = interp_nuclide_xs(name, temp, bg_xs_1);
  auto xs_2 = interp_nuclide_xs(name, temp, bg_xs_2);

  xt::xtensor<double, 1> Etr = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 2> Es_tr = xt::zeros<double>({ngroups_, ngroups_});
  xt::xtensor<double, 1> Ef = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> vEf = xt::zeros<double>({ngroups_});
  xt::xtensor<double, 1> chi = xt::zeros<double>({ngroups_});

  double vEf_sum_1 = 0.;
  double vEf_sum_2 = 0.;
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
    Ea(g) = f1_g * xs_1->Ea(g) + f2_g * xs_2->Ea(g);
    Ef(g) = f1_g * xs_1->Ef(g) + f2_g * xs_2->Ef(g);
    for (std::size_t g_out = 0; g_out < ngroups_; g_out++) {
      Es_tr(g, g_out) =
          f1_g * xs_1->Es_tr(g, g_out) + f2_g * xs_2->Es_tr(g, g_out);
    }
    Etr(g) = Ea(g) + xt::sum(xt::view(Es_tr, g, xt::all()))();

    const double vEf1 = f1_g * xs_1->vEf(g);
    const double vEf2 = f2_g * xs_2->vEf(g);
    vEf(g) = vEf1 + vEf2;
    vEf_sum_1 += vEf1;
    vEf_sum_2 += vEf2;
  }

  if (vEf_sum_1 + vEf_sum_2 > 0.) {
    double chi_sum = 0;
    for (std::size_t g = 0; g < ngroups_; g++) {
      chi(g) = (vEf_sum_1 * xs_1->chi(g) + vEf_sum_2 * xs_2->chi(g)) /
               (vEf_sum_1 + vEf_sum_2);
      chi_sum += chi(g);
    }
    if (chi_sum > 0.) chi /= chi_sum;
  }

  std::shared_ptr<CrossSection> xs_out =
      std::make_shared<CrossSection>(Etr, Ea, Es_tr, Ef, vEf, chi);
  *xs_out *= N;
  return xs_out;
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
                          const xt::xtensor<double, 2> nE, std::size_t it,
                          double f_temp) const {
  if (f_temp > 0.) {
    E = (1. - f_temp) * xt::view(nE, it, xt::all()) +
        f_temp * xt::view(nE, it + 1, xt::all());
  } else {
    E = xt::view(nE, it, xt::all());
  }
}

void NDLibrary::interp_1d(xt::xtensor<double, 1>& E,
                          const xt::xtensor<double, 3> nE, std::size_t it,
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
                          const xt::xtensor<double, 4> nE, std::size_t it,
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

}  // namespace scarabee

// References
// [1] H. Koike, K. Yamaji, K. Kirimura, D. Sato, H. Matsumoto, and A. Yamamoto,
//     “Advanced resonance self-shielding method for gray resonance treatment in
//     lattice physics code GALAXY,” J. Nucl. Sci. Technol., vol. 49, no. 7,
//     pp. 725–747, 2012, doi: 10.1080/00223131.2012.693885.
//
// [2] N. Gibson, “Novel Resonance Self-Shielding Methods for Nuclear Reactor
//     Analysis,” Massachusetts Institute of Technology, 2016.