#ifndef CELL_H
#define CELL_H

#include <cross_section.hpp>
#include <moc/direction.hpp>
#include <moc/segment.hpp>
#include <moc/vector.hpp>
#include <moc/flat_source_region.hpp>
#include <moc/surface.hpp>
#include <utils/constants.hpp>

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace scarabee {

class Cell {
 public:
  bool inside(const Vector& r, const Direction& u) const;

  double distance(const Vector& r, const Direction& u) const;

  UniqueFSR get_fsr(const Vector& r, const Direction& u) const;

  std::size_t num_fsrs() const { return fsrs_.size(); }

  std::size_t ngroups() const { return fsrs_.front().xs()->ngroups(); }

  std::set<std::size_t> get_all_fsr_ids() const;

  std::size_t get_num_fsr_instances(std::size_t id) const;

  void fill_fsrs(std::map<std::size_t, const FlatSourceRegion*>& fsrs) const;

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
