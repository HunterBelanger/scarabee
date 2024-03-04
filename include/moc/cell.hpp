#ifndef CELL_H
#define CELL_H

#include <transport_xs.hpp>
#include <moc/direction.hpp>
#include <moc/segment.hpp>
#include <moc/vector.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/surface.hpp>
#include <utils/constants.hpp>

#include <memory>
#include <vector>

class Cell {
 public:
  Cell(std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
       std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax);

  std::vector<Segment> trace_segments(Vector& r, const Direction& u);

  bool inside(const Vector& r, const Direction& u) const;
  double distance(const Vector& r, const Direction& u) const;

  FlatSourceRegion& get_fsr(const Vector& r, const Direction& u);
  const FlatSourceRegion& get_fsr(const Vector& r, const Direction& u) const;

 protected:
  std::vector<FlatSourceRegion> fsrs_;
  std::shared_ptr<Surface> x_min_, y_min_, x_max_, y_max_;
};

#endif
