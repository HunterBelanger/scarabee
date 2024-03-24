#ifndef MOC_DRIVER_H
#define MOC_DRIVER_H

#include <moc/cartesian_2d.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/track.hpp>

#include <memory>
#include <vector>

enum class BoundaryCondition : std::uint8_t { Reflective, Vacuum };

class MOCDriver {
 public:
  MOCDriver(std::shared_ptr<Cartesian2D> geometry,
            BoundaryCondition xmin = BoundaryCondition::Reflective,
            BoundaryCondition xmax = BoundaryCondition::Reflective,
            BoundaryCondition ymin = BoundaryCondition::Reflective,
            BoundaryCondition ymax = BoundaryCondition::Reflective);

  bool drawn() const { return !angle_info_.empty(); }

  void draw_tracks(std::uint32_t n_angles, double d);

  BoundaryCondition& x_min_bc() { return x_min_bc_; }
  const BoundaryCondition& x_min_bc() const { return x_min_bc_; }

  BoundaryCondition& x_max_bc() { return x_max_bc_; }
  const BoundaryCondition& x_max_bc() const { return x_max_bc_; }

  BoundaryCondition& y_min_bc() { return y_min_bc_; }
  const BoundaryCondition& y_min_bc() const { return y_min_bc_; }

  BoundaryCondition& y_max_bc() { return y_max_bc_; }
  const BoundaryCondition& y_max_bc() const { return y_max_bc_; }

 private:
  struct AngleInfo {
    double phi;        // Azimuthal angle for track
    double d;          // Spacing for trackings of this angle
    double wgt;        // Weight for tracks with this angle
    std::uint32_t nx;  // Number of tracks starting on the -y boundary
    std::uint32_t ny;  // Number of tracks starting on the -x boundary
  };

  std::vector<AngleInfo> angle_info_;       // Information for all angles
  std::vector<std::vector<Track>> tracks_;  // All tracks, indexed by angle
  std::vector<FlatSourceRegion*> fsrs_;     // All FSRs in the geometry
  std::shared_ptr<Cartesian2D> geometry_;   // Geometry for the problem
  BoundaryCondition x_min_bc_, x_max_bc_, y_min_bc_, y_max_bc_;
};

#endif