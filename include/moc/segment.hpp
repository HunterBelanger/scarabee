#ifndef SEGMENT_H
#define SEGMENT_H

#include <transport_xs.hpp>
#include <moc/flat_source_region.hpp>

#include <xtensor/xarray.hpp>

class Segment {
 public:
  Segment(FlatSourceRegion* fsr, double length) : fsr_(fsr), length_(length) {}

  double length() const { return length_; }

  const TransportXS& xs() const { return *fsr_->xs(); }

  xt::xarray<double>& flux() { return fsr_->flux(); }
  const xt::xarray<double>& flux() const { return fsr_->flux(); }

  xt::xarray<double>& source() { return fsr_->source(); }
  const xt::xarray<double>& source() const { return fsr_->source(); }

 private:
  FlatSourceRegion* fsr_;
  double length_;
};

#endif
