#ifndef SCARABEE_CRITICALITY_SPECTRUM_H
#define SCARABEE_CRITICALITY_SPECTRUM_H

#include <cross_section.hpp>

#include <xtensor/xtensor.hpp>

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
};

class P1CriticalitySpectrum : public CriticalitySpectrum {
 public:
  P1CriticalitySpectrum(std::shared_ptr<CrossSection> xs);
};

class B1CriticalitySpectrum : public CriticalitySpectrum {
 public:
  B1CriticalitySpectrum(std::shared_ptr<CrossSection> xs);
};

}  // namespace scarabee

#endif