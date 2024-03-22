#ifndef SEGMENT_H
#define SEGMENT_H

#include <transport_xs.hpp>
#include <moc/flat_source_region.hpp>

#include <vector>

class Segment {
 public:
  Segment(FlatSourceRegion* fsr, double length) : fsr_(fsr), length_(length) {}

  double length() const { return length_; }

  const TransportXS& xs() const { return *fsr_->xs(); }

  std::vector<double>& flux() { return fsr_->flux(); }
  const std::vector<double>& flux() const { return fsr_->flux(); }

  std::vector<double>& source() { return fsr_->source(); }
  const std::vector<double>& source() const { return fsr_->source(); }

 private:
  FlatSourceRegion* fsr_;
  double length_;
};

#endif
