#ifndef SCARABEE_CYLINDRICAL_FLUX_SOLVER_H
#define SCARABEE_CYLINDRICAL_FLUX_SOLVER_H

#include <cylindrical_cell.hpp>

#include <ndarray.hpp>

#include <memory>
#include <vector>

class CylindricalFluxSolver {
 public:
  CylindricalFluxSolver(std::shared_ptr<CylindricalCell> cell);

  std::uint32_t ngroups() const { return cell_->ngroups(); }
  std::size_t nregions() const { return cell_->nregions(); }

  double flux(std::uint32_t g, std::size_t i) const { return flux_(g, i); }

  double Q(std::uint32_t g, std::size_t i) const;

  double keff() const { return k_; }
  double keff_tolerance() const { return k_tol_; }
  void set_keff_tolerance(double ktol);

  double albedo() const { return a_; }
  void set_albedo(double a);

  double j_ext(std::uint32_t g) const { return j_ext_[g]; }
  void set_j_ext(std::uint32_t g, double j) { j_ext_[g] = j; }

  void solve();
  bool solved() const { return solved_; }

 private:
  NDArray<double> flux_;
  std::vector<double> j_ext_;
  std::shared_ptr<CylindricalCell> cell_;
  double k_;
  double a_;
  double k_tol_;
  bool solved_;

  double calc_keff(const NDArray<double>& flux) const;
};

#endif
