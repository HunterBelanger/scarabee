#include <moc/track.hpp>
#include <utils/scarabee_exception.hpp>

Track::Track(const Vector& r_start, const Vector& r_end, const Direction& u,
             double phi, double wgt, const std::vector<Segment>& segments)
    : segments_(segments),
      r_start_(r_start),
      r_end_(r_end),
      u_(u),
      entry_track_(nullptr),
      exit_track_(nullptr),
      weight_(wgt),
      phi_(phi) {}

Segment& Track::at(std::size_t i) {
  if (i >= this->size()) {
    throw ScarabeeException("Index out of range.");
  }

  return (*this)[i];
}

const Segment& Track::at(std::size_t i) const {
  if (i >= this->size()) {
    throw ScarabeeException("Index out of range.");
  }

  return (*this)[i];
}