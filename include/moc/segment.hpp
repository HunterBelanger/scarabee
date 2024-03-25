#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <transport_xs.hpp>

#include <xtensor/xarray.hpp>

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

 private:
  FlatSourceRegion* fsr_;
  double length_;
};

#endif
