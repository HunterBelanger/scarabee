#include <moc/track.hpp>
#include <utils/scarabee_exception.hpp>

Track::Track(const Vector& start, const Vector& end, const Direction& dir,
             double phi, double wgt, const std::vector<Segment>& segments)
    : entry_flux_(),
      exit_flux_(),
      segments_(segments),
      start_(start),
      end_(end),
      dir_(dir),
      entry_track_(nullptr),
      exit_track_(nullptr),
      weight_(wgt),
      phi_(phi),
      start_bc_(BoundaryCondition::Reflective),
      end_bc_(BoundaryCondition::Reflective) {}

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