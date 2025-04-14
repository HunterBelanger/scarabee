#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <data/cross_section.hpp>

#include <xtensor/containers/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>

#include <memory>
#include <optional>

namespace scarabee {

class Segment {
 public:
  Segment(const FlatSourceRegion* fsr, double length, std::size_t indx)
      : xs_(fsr->xs()),
        volume_(fsr->volume()),
        length_(length),
        fsr_indx_(indx) {}

  // Here for use with cereal and std::vector
  Segment() {}

  double length() const { return length_; }

  void set_length(double l) { length_ = l; }

  double volume() const { return volume_; }

  const std::shared_ptr<CrossSection>& xs() const { return xs_; }

  std::size_t fsr_indx() const { return fsr_indx_; }

 private:
  std::shared_ptr<CrossSection> xs_;
  double volume_;
  double length_;
  std::size_t fsr_indx_;

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(xs_), CEREAL_NVP(volume_), CEREAL_NVP(length_),
        CEREAL_NVP(fsr_indx_));
  }
};

}  // namespace scarabee

#endif
