#ifndef CELL_H
#define CELL_H

#include <transport_xs.hpp>
#include <moc/direction.hpp>
#include <moc/vector.hpp>

#include <xtensor/xarray.hpp>

#include <memory>

class Cell {
 public:
  Cell(double xmin, double xmax, double ymin, double ymax);
  virtual ~Cell() = default;

  virtual void trace_segments(const Vector& r, const Direction& u) const = 0;
  virtual void get_fsr(const Vector& r, const Direction& u) const = 0;
  virtual const TransportXS* get_xs(const Vector& r, const Direction& u) const = 0;

  double x_min() const { return x_min_; }
  double x_max() const { return x_max_; }
  double y_min() const { return y_min_; }
  double y_max() const { return y_max_; }

  void set_x_min(const double& val) {
    x_min_ = val;
    check_values();
  }
  void set_x_max(const double& val) {
    x_max_ = val;
    check_values();
  }
  void set_y_min(const double& val) {
    y_min_ = val;
    check_values();
  }
  void set_y_max(const double& val) {
    y_max_ = val;
    check_values();
  }

 private:
  xt::xarray<double> fsrs_;
  xt::xarray<std::shared_ptr<TransportXS>> xss_;
  double x_min_, x_max_;
  double y_min_, y_max_;

  virtual void check_values() const;
};

#endif
