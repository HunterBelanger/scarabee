#ifndef TRACK_H
#define TRACK_H

#include <moc/segment.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>
#include <moc/boundary_condition.hpp>
#include <utils/constants.hpp>

#include <xtensor/xtensor.hpp>

#include <vector>

class Track {
 public:
  Track(const Vector& start, const Vector& end, const Direction& dir,
        double phi, double wgt, const std::vector<Segment>& segments);

  double weight() const { return weight_; }
  double phi() const { return phi_; }
  const std::vector<Segment>& segments() const { return segments_; }

  const Vector& start_pos() const { return start_; }
  const Vector& end_pos() const { return end_; }
  const Direction& dir() const { return dir_; } 

  BoundaryCondition& start_bc() { return start_bc_; }
  const BoundaryCondition& start_bc() const { return start_bc_; }

  BoundaryCondition& end_bc() { return end_bc_; }
  const BoundaryCondition& end_bc() const { return end_bc_; }

  xt::xtensor<double, 2>& entry_flux() { return entry_flux_; }
  const xt::xtensor<double, 2>& entry_flux() const { return entry_flux_; }

  xt::xtensor<double, 2>& exit_flux() { return exit_flux_; }
  const xt::xtensor<double, 2>& exit_flux() const { return exit_flux_; }

  Track* entry_track() { return entry_track_; }
  void set_entry_track(Track* t) { entry_track_ = t; }

  Track* exit_track() { return exit_track_; }
  void set_exit_track(Track* t) { exit_track_ = t; }

  // Values and methods for the reverse direction
  const Vector& rstart_pos() const { return end_; }
  const Vector& rend_pos() const { return start_; }
  Direction rdir() const { return -dir_; }
  double rphi() const { return PI + phi_; } 

  xt::xtensor<double, 2>& rentry_flux() { return exit_flux_; }
  const xt::xtensor<double, 2>& rentry_flux() const { return exit_flux_; }

  xt::xtensor<double, 2>& rexit_flux() { return entry_flux_; }
  const xt::xtensor<double, 2>& rexit_flux() const { return entry_flux_; }

  Track* rentry_track() { return exit_track_; }
  void set_rentry_track(Track* t) { exit_track_ = t; }

  Track* rexit_track() { return entry_track_; }
  void set_rexit_track(Track* t) { entry_track_ = t; }

  // Indexing is only done in forward direction
  std::size_t size() const { return segments_.size(); }

  Segment& operator[](std::size_t i) { return segments_[i]; }
  const Segment& operator[](std::size_t i) const { return segments_[i]; }

  Segment& at(std::size_t i);
  const Segment& at(std::size_t i) const;

  //--------------------------------------------------------------------------
  // Iterators
  using iterator = std::vector<Segment>::iterator;
  using const_iterator = std::vector<Segment>::const_iterator;
  using reverse_iterator = std::vector<Segment>::reverse_iterator;
  using const_reverse_iterator = std::vector<Segment>::const_reverse_iterator;

  iterator begin() { return segments_.begin(); }
  const_iterator begin() const { return segments_.begin(); }
  const_iterator cbegin() const { return segments_.cbegin(); }

  iterator end() { return segments_.end(); }
  const_iterator end() const { return segments_.end(); }
  const_iterator cend() const { return segments_.cend(); }

  reverse_iterator rbegin() { return segments_.rbegin(); }
  const_reverse_iterator rbegin() const { return segments_.rbegin(); }
  const_reverse_iterator crbegin() const { return segments_.crbegin(); }

  reverse_iterator rend() { return segments_.rend(); }
  const_reverse_iterator rend() const { return segments_.rend(); }
  const_reverse_iterator crend() const { return segments_.crend(); }

 private:
  xt::xtensor<double, 2> entry_flux_; // Indexed on group and polar angle
  xt::xtensor<double, 2> exit_flux_;
  std::vector<Segment> segments_;
  Vector start_;
  Vector end_;
  Direction dir_;
  Track* entry_track_;
  Track* exit_track_;
  double weight_;  // Track weight
  double phi_;     // Azimuthal angle of Track
  BoundaryCondition start_bc_;
  BoundaryCondition end_bc_;
};

#endif