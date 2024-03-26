#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <transport_xs.hpp>

#include <xtensor/xtensor.hpp>

#include <vector>

class Segment {
 public:
  Segment(FlatSourceRegion* fsr, double length) : fsr_(fsr), length_(length) {}

  double length() const { return length_; }

  const TransportXS& xs() const { return *fsr_->xs(); }

  xt::xtensor<double, 1>& flux() { return fsr_->flux(); }
  const xt::xtensor<double, 1>& flux() const { return fsr_->flux(); }

  xt::xtensor<double, 1>& source() { return fsr_->source(); }
  const xt::xtensor<double, 1>& source() const { return fsr_->source(); }

  xt::xtensor<double, 2>& exp() { return exp_; }
  const xt::xtensor<double, 2>& exp() const { return exp_; }

 private:
  xt::xtensor<double, 2> exp_; // exp (- Etg * l / sin(theta))
  FlatSourceRegion* fsr_;
  double length_;
};

#endif
