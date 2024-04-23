#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <transport_xs.hpp>

#include <xtensor/xtensor.hpp>

#include <memory>

namespace scarabee {

class Segment {
 public:
  Segment(const FlatSourceRegion* fsr, double length, std::size_t indx)
      : fsr_(fsr), length_(length), fsr_indx_(indx) {}

  double length() const { return length_; }

  void set_length(double l) { length_ = l; }

  double volume() const { return fsr_->volume(); }

  const std::shared_ptr<TransportXS>& xs() const { return fsr_->xs(); }

  std::size_t fsr_indx() const { return fsr_indx_; }

 private:
  const FlatSourceRegion* fsr_;
  double length_;
  std::size_t fsr_indx_;
};

}  // namespace scarabee

#endif
