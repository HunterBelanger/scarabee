#include <diffusion/diffusion_data.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xmanipulation.hpp>
#include <xtensor-io/xnpz.hpp>

#include <filesystem>
#include <memory>

namespace scarabee {

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                             const xt::xtensor<double, 2>& form_factors)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  this->set_form_factors(form_factors);
}

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                             const xt::xtensor<double, 2>& form_factors,
                             const xt::xtensor<double, 2>& adf,
                             const xt::xtensor<double, 2>& cdf)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  this->set_form_factors(form_factors);
  this->set_adf(adf);
  this->set_cdf(cdf);
}

void DiffusionData::set_form_factors(const xt::xtensor<double, 2>& ff) {
  for (std::size_t i = 0; i < ff.size(); i++) {
    if (ff.flat(i) < 0.) {
      auto mssg = "Form factors must be positive.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  form_factors_ = ff;
}

void DiffusionData::set_adf(const xt::xtensor<double, 2>& adf) {
  if (adf.shape()[0] != this->ngroups()) {
    auto mssg =
        "Number of groups in ADF array does not agree with the number of "
        "groups in the diffusion cross section.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (adf.shape()[1] != static_cast<std::size_t>(4)) {
    auto mssg = "The ADF array must have four columns.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < adf.size(); i++) {
    if (adf.flat(i) < 0.) {
      auto mssg = "ADF values must be positive.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  adf_ = adf;
}

void DiffusionData::set_cdf(const xt::xtensor<double, 2>& cdf) {
  if (cdf.shape()[0] != this->ngroups()) {
    auto mssg =
        "Number of groups in CDF array does not agree with the number of "
        "groups in the diffusion cross section.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (cdf.shape()[1] != static_cast<std::size_t>(4)) {
    auto mssg = "The CDF array must have four columns.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < cdf.size(); i++) {
    if (cdf.flat(i) < 0.) {
      auto mssg = "CDF values must be positive.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  cdf_ = cdf;
}

void DiffusionData::rotate_clockwise() {
  if (form_factors_.size() > 0) {
    xt::xtensor<double, 2> temp;

    // First transpose the matrix
    temp = xt::transpose(form_factors_);

    // Now we reverse each row
    form_factors_ = xt::flip(temp, 1);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = adf_(g, ADF::XN);

      adf_(g, ADF::XN) = adf_(g, ADF::YN);
      adf_(g, ADF::YN) = adf_(g, ADF::XP);
      adf_(g, ADF::XP) = adf_(g, ADF::YP);
      adf_(g, ADF::YP) = temp2;
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = cdf_(g, CDF::I);

      cdf_(g, CDF::I) = cdf_(g, CDF::II);
      cdf_(g, CDF::II) = cdf_(g, CDF::III);
      cdf_(g, CDF::III) = cdf_(g, CDF::IV);
      cdf_(g, CDF::IV) = temp2;
    }
  }
}

void DiffusionData::rotate_counterclockwise() {
  if (form_factors_.size() > 0) {
    xt::xtensor<double, 2> temp;

    // First we reverse each row
    temp = xt::flip(form_factors_, 1);

    // Now transpose the matrix
    form_factors_ = xt::transpose(temp);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = adf_(g, ADF::XN);

      adf_(g, ADF::XN) = adf_(g, ADF::YP);
      adf_(g, ADF::YP) = adf_(g, ADF::XP);
      adf_(g, ADF::XP) = adf_(g, ADF::YN);
      adf_(g, ADF::YP) = temp2;
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = cdf_(g, CDF::I);

      cdf_(g, CDF::I) = cdf_(g, CDF::IV);
      cdf_(g, CDF::IV) = cdf_(g, CDF::III);
      cdf_(g, CDF::III) = cdf_(g, CDF::II);
      cdf_(g, CDF::II) = temp2;
    }
  }
}

void DiffusionData::save(const std::string& fname) const {
  if (std::filesystem::exists(fname)) {
    std::filesystem::remove(fname);
  }

  const std::size_t NG = this->ngroups();

  xt::xtensor<double, 2> xs_data = xt::zeros<double>({5 + NG, NG});
  for (std::size_t g = 0; g < NG; g++) {
    xs_data(0, g) = this->D(g);
    xs_data(1, g) = this->Ea(g);
    xs_data(2, g) = this->Ef(g);
    xs_data(3, g) = this->vEf(g);
    xs_data(4, g) = this->chi(g);

    // Save Es
    for (std::size_t gout = 0; gout < NG; gout++) {
      xs_data(5 + g, gout) = this->Es(g, gout);
    }
  }

  // Save the xs data to the npz
  xt::dump_npz(fname, "xs", xs_data);

  if (form_factors_.size() > 0) {
    xt::dump_npz(fname, "form_factors", form_factors_);
  }

  if (adf_.size() > 0) {
    xt::dump_npz(fname, "adf", adf_);
  }

  if (cdf_.size() > 0) {
    xt::dump_npz(fname, "cdf", cdf_);
  }
}

std::shared_ptr<DiffusionData> DiffusionData::load(const std::string& fname) {
  auto npz = xt::load_npz(fname);

  // First, read in and create the DiffsuionCrossSection
  if (npz.find("xs") == npz.end()) {
    auto mssg = "No \"xs\" entry in the DiffusionData NPZ \"" + fname + "\".";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  xt::xtensor<double, 2> xs_data = npz.at("xs").cast<double>();
  const std::size_t NG = xs_data.shape()[1];
  xt::xtensor<double, 1> D = xt::view(xs_data, 0, xt::all());
  xt::xtensor<double, 1> Ea = xt::view(xs_data, 1, xt::all());
  xt::xtensor<double, 1> Ef = xt::view(xs_data, 2, xt::all());
  xt::xtensor<double, 1> vEf = xt::view(xs_data, 3, xt::all());
  xt::xtensor<double, 1> chi = xt::view(xs_data, 4, xt::all());
  xt::xtensor<double, 2> Es =
      xt::view(xs_data, xt::range(5, 5 + NG), xt::all());
  auto xs = std::make_shared<DiffusionCrossSection>(D, Ea, Es, Ef, vEf, chi);

  std::unique_ptr<xt::xtensor<double, 2>> form_factors, adf, cdf;

  // Load form factors
  if (npz.find("form_factors") != npz.end()) {
    form_factors = std::make_unique<xt::xtensor<double, 2>>();
    *form_factors = npz.at("form_factors").cast<double>();
  }

  // Load ADF
  if (npz.find("adf") != npz.end()) {
    adf = std::make_unique<xt::xtensor<double, 2>>();
    *adf = npz.at("adf").cast<double>();
  }

  // Load CDF
  if (npz.find("cdf") != npz.end()) {
    cdf = std::make_unique<xt::xtensor<double, 2>>();
    *cdf = npz.at("cdf").cast<double>();
  }

  // Return data
  if (form_factors && adf && cdf) {
    return std::make_shared<DiffusionData>(xs, *form_factors, *adf, *cdf);
  } else if (form_factors && adf) {
    auto out = std::make_shared<DiffusionData>(xs, *form_factors);
    out->set_adf(*adf);
    return out;
  } else if (form_factors && cdf) {
    auto out = std::make_shared<DiffusionData>(xs, *form_factors);
    out->set_cdf(*cdf);
    return out;
  } else if (form_factors) {
    return std::make_shared<DiffusionData>(xs, *form_factors);
  } else if (adf && cdf) {
    auto out = std::make_shared<DiffusionData>(xs);
    out->set_adf(*adf);
    out->set_cdf(*cdf);
    return out;
  } else if (adf) {
    auto out = std::make_shared<DiffusionData>(xs);
    out->set_adf(*adf);
    return out;
  } else if (cdf) {
    auto out = std::make_shared<DiffusionData>(xs);
    out->set_cdf(*cdf);
    return out;
  } else {
    return std::make_shared<DiffusionData>(xs);
  }

  // SHOULD NEVER GET HERE
  auto mssg = "Something went wrong loading DiffusionData.";
  spdlog::error(mssg);
  throw ScarabeeException(mssg);
  return nullptr;
}

}  // namespace scarabee
