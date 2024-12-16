#include <diffusion/diffusion_data.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <filesystem>
#include <fstream>
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

  std::ofstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryOutputArchive arc(file);

  arc(*this);
}

std::shared_ptr<DiffusionData> DiffusionData::load(const std::string& fname) {
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  std::shared_ptr<DiffusionData> out(new DiffusionData());

  std::ifstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryInputArchive arc(file);

  arc(*out);

  return out;
}

}  // namespace scarabee
