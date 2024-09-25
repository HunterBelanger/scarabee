#include <moc/track.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

namespace scarabee {

Track::Track(const Vector& entry, const Vector& exit, const Direction& dir,
             double phi, double wgt, double width,
             const std::vector<Segment>& segments,
             std::size_t forward_phi_index, std::size_t backward_phi_index)
    : entry_flux_(),
      exit_flux_(),
      segments_(segments),
      entry_(entry),
      exit_(exit),
      dir_(dir),
      entry_track_flux_(nullptr),
      exit_track_flux_(nullptr),
      wgt_(wgt),
      width_(width),
      phi_(phi),
      entry_bc_(BoundaryCondition::Reflective),
      exit_bc_(BoundaryCondition::Reflective),
      forward_phi_index_(forward_phi_index),
      backward_phi_index_(backward_phi_index) {}

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
