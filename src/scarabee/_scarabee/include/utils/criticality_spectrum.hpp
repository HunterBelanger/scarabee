#ifndef SCARABEE_CRITICALITY_SPECTRUM_H
#define SCARABEE_CRITICALITY_SPECTRUM_H

#include <data/cross_section.hpp>
#include <utils/serialization.hpp>

#include <xtensor/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>

namespace scarabee {

class CriticalitySpectrum {
 public:
  std::size_t ngroups() const { return flux_.size(); }

  double k_inf() const { return k_inf_; }

  double B2() const { return B2_; }
  double buckling() const { return B2(); }

  const xt::xtensor<double, 1>& flux() const { return flux_; }
  const xt::xtensor<double, 1>& current() const { return current_; }
  const xt::xtensor<double, 1>& diff_coeff() const { return diff_coeff_; }

  double flux(std::size_t g) const { return flux_(g); }
  double current(std::size_t g) const { return current_(g); }
  double diff_coeff(std::size_t g) const { return diff_coeff_(g); }

 protected:
  double k_inf_, B2_;
  xt::xtensor<double, 1> flux_;
  xt::xtensor<double, 1> current_;
  xt::xtensor<double, 1> diff_coeff_;

  CriticalitySpectrum() : k_inf_(), B2_(), flux_(), current_(), diff_coeff_() {}

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(k_inf_), CEREAL_NVP(B2_), CEREAL_NVP(flux_),
        CEREAL_NVP(current_), CEREAL_NVP(diff_coeff_));
  }
};

class P1CriticalitySpectrum : public CriticalitySpectrum {
 public:
  P1CriticalitySpectrum(std::shared_ptr<CrossSection> xs);
  P1CriticalitySpectrum(std::shared_ptr<CrossSection> xs, double B);

 private:
  friend class cereal::access;
  P1CriticalitySpectrum() {}
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<CriticalitySpectrum>(this));
  }
};

class B1CriticalitySpectrum : public CriticalitySpectrum {
 public:
  B1CriticalitySpectrum(std::shared_ptr<CrossSection> xs);
  B1CriticalitySpectrum(std::shared_ptr<CrossSection> xs, double B);

 private:
  friend class cereal::access;
  B1CriticalitySpectrum() {}
  template <class Archive>
  void serialize(Archive& arc) {
    arc(cereal::base_class<CriticalitySpectrum>(this));
  }
};

}  // namespace scarabee

#endif