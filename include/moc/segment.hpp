#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <transport_xs.hpp>

#include <xtensor/xtensor.hpp>

namespace scarabee {

class Segment {
 public:
  Segment(FlatSourceRegion* fsr, double length) : fsr_(fsr), length_(length) {}

  double length() const { return length_; }

  void set_length(double l) { length_ = l; }

  double volume() const { return fsr_->volume(); }

  const TransportXS& xs() const { return *fsr_->xs(); }

  std::size_t fsr_indx() const { return fsr_->indx(); }

  xt::xtensor<double, 2>& exp() { return exp_; }
  const xt::xtensor<double, 2>& exp() const { return exp_; }

 private:
  xt::xtensor<double, 2> exp_;  // exp (- Etg * l / sin(theta))
  FlatSourceRegion* fsr_;
  double length_;
};

}  // namespace scarabee

#endif
