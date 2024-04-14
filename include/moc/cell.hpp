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

namespace scarabee {

class Cell {
 public:
  virtual ~Cell() = default;

  virtual std::shared_ptr<Cell> clone() const = 0;

  std::vector<Segment> trace_segments(Vector& r, const Direction& u);
  double trace_segments(Vector& r, const Direction& u, std::vector<Segment>& segments);

  bool inside(const Vector& r, const Direction& u) const;

  double distance(const Vector& r, const Direction& u) const;

  FlatSourceRegion& get_fsr(const Vector& r, const Direction& u);
  const FlatSourceRegion& get_fsr(const Vector& r, const Direction& u) const;

  std::size_t num_fsrs() const { return fsrs_.size(); }

  void append_fsrs(std::vector<FlatSourceRegion*>& fsrs) {
    for (auto& fsr : fsrs_) fsrs.push_back(&fsr);
  }

  double dx() const { return x_max_->x0() - x_min_->x0(); }
  double dy() const { return y_max_->y0() - y_min_->y0(); }

 protected:
  std::vector<FlatSourceRegion> fsrs_;
  std::shared_ptr<Surface> x_min_, y_min_, x_max_, y_max_;

  Cell(double dx, double dy);
  void check_surfaces() const;
};

}  // namespace scarabee

#endif
