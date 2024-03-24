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
  std::vector<Segment> trace_segments(Vector& r, const Direction& u);
  void trace_segments(Vector& r, const Direction& u,
                      std::vector<Segment>& segments);

  bool inside(const Vector& r, const Direction& u) const;

  double distance(const Vector& r, const Direction& u) const;

  FlatSourceRegion& get_fsr(const Vector& r, const Direction& u);
  const FlatSourceRegion& get_fsr(const Vector& r, const Direction& u) const;

  std::size_t num_fsrs() const { return fsrs_.size(); }

  void append_fsrs(std::vector<FlatSourceRegion*>& fsrs) {
    for (auto& fsr : fsrs_) fsrs.push_back(&fsr);
  }

  const std::shared_ptr<Surface>& x_min() const { return x_min_; }
  const std::shared_ptr<Surface>& x_max() const { return x_max_; }
  const std::shared_ptr<Surface>& y_min() const { return y_min_; }
  const std::shared_ptr<Surface>& y_max() const { return y_max_; }

  void set_x_min(const std::shared_ptr<Surface>& xm);
  void set_x_max(const std::shared_ptr<Surface>& xm);
  void set_y_min(const std::shared_ptr<Surface>& ym);
  void set_y_max(const std::shared_ptr<Surface>& ym);

 protected:
  std::vector<FlatSourceRegion> fsrs_;
  std::shared_ptr<Surface> x_min_, y_min_, x_max_, y_max_;

  Cell(std::shared_ptr<Surface>& xmin, std::shared_ptr<Surface>& xmax,
       std::shared_ptr<Surface>& ymin, std::shared_ptr<Surface>& ymax);
  void check_surfaces() const;
};

#endif
