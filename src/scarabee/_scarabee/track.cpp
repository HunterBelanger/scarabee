#include <moc/track.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

namespace scarabee {

Track::Track(const Vector& entry, const Vector& exit, const Direction& dir,
             double phi, double wgt, const std::vector<Segment>& segments)
    : entry_flux_(),
      exit_flux_(),
      segments_(segments),
      entry_(entry),
      exit_(exit),
      dir_(dir),
      entry_track_flux_(nullptr),
      exit_track_flux_(nullptr),
      weight_(wgt),
      phi_(phi),
      entry_bc_(BoundaryCondition::Reflective),
      exit_bc_(BoundaryCondition::Reflective) {}

Segment& Track::at(std::size_t i) {
  if (i >= this->size()) {
    auto mssg = "Index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return (*this)[i];
}

const Segment& Track::at(std::size_t i) const {
  if (i >= this->size()) {
    auto mssg = "Index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return (*this)[i];
}

}  // namespace scarabee
