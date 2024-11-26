#ifndef SEGMENT_H
#define SEGMENT_H

#include <moc/flat_source_region.hpp>
#include <data/cross_section.hpp>

#include <xtensor/xtensor.hpp>

#include <memory>
#include <optional>

namespace scarabee {

class Segment {
 public:
  Segment(const FlatSourceRegion* fsr, double length, std::size_t indx)
      : fsr_(fsr),
        length_(length),
        fsr_indx_(indx),
        entry_cmfd_surface_(std::nullopt),
        exit_cmfd_surface_(std::nullopt) {}

  double length() const { return length_; }

  void set_length(double l) { length_ = l; }

  double volume() const { return fsr_->volume(); }

  const std::shared_ptr<CrossSection>& xs() const { return fsr_->xs(); }

  std::size_t fsr_indx() const { return fsr_indx_; }

  std::optional<std::size_t>& entry_cmfd_surface() {
    return entry_cmfd_surface_;
  }
  const std::optional<std::size_t>& entry_cmfd_surface() const {
    return entry_cmfd_surface_;
  }

  std::optional<std::size_t>& exit_cmfd_surface() { return exit_cmfd_surface_; }
  const std::optional<std::size_t>& exit_cmfd_surface() const {
    return exit_cmfd_surface_;
  }

 private:
  const FlatSourceRegion* fsr_;
  double length_;
  std::size_t fsr_indx_;
  std::optional<std::size_t> entry_cmfd_surface_;
  std::optional<std::size_t> exit_cmfd_surface_;
};

}  // namespace scarabee

#endif
