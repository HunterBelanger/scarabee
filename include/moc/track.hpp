#ifndef TRACK_H
#define TRACK_H

#include <moc/segment.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>

#include <vector>

class Track {
 public:
  Track(const Vector& r_start, const Vector& r_end, const Direction& u,
        double phi, double wgt, const std::vector<Segment>& segments);

  double weight() const { return weight_; }
  double phi() const { return phi_; }
  const std::vector<Segment>& segments() const { return segments_; }
  const Vector& r_start() const { return r_start_; }
  const Vector& r_end() const { return r_end_; }
  const Direction& u() const { return u_; }

  Track* entry_track() { return entry_track_; }
  void set_entry_track(Track* t) { entry_track_ = t; }

  Track* exit_track() { return exit_track_; }
  void set_exit_track(Track* t) { exit_track_ = t; }

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
  std::vector<Segment> segments_;
  Vector r_start_;
  Vector r_end_;
  Direction u_;
  Track* entry_track_;
  Track* exit_track_;
  double weight_;  // Track weight
  double phi_;     // Azimuthal angle of Track
};

#endif