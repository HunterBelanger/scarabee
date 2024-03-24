#ifndef MOC_DRIVER_H
#define MOC_DRIVER_H

#include <moc/cartesian_2d.hpp>
#include <moc/track.hpp>

#include <memory>
#include <vector>

class MOCDriver {
 public:
  MOCDriver(std::shared_ptr<Cartesian2D> geometry, double d, std::uint32_t na);

 private:
  struct AngleInfo {
    double phi;        // Azimuthal angle for track
    double d;          // Spacing for trackings of this angle
    double wgt;        // Weight for tracks with this angle
    std::uint32_t nx;  // Number of tracks starting on the -y boundary
    std::uint32_t ny;  // Number of tracks starting on the -x boundary
  };

  std::vector<AngleInfo> angle_info_;  // Information for all angles
  std::vector<std::vector<Track>> tracks_;
  std::shared_ptr<Cartesian2D> geometry_;
  std::uint32_t n_angles_;  // Number of azimuthal angles
};

#endif