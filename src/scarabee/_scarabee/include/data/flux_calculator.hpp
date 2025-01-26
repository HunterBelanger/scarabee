#ifndef SCARABEE_FLUX_CALCULATOR_H
#define SCARABEE_FLUX_CALCULATOR_H

#include <xtensor/xtensor.hpp>

#include <cstdint>
#include <vector>

namespace scarabee {

class FluxCalculator {
 public:
  struct BackgroundNuclide {
    double A;
    double a;
    double sig_b;
    std::size_t g_max;
  };

  FluxCalculator(const xt::xtensor<double, 1>& energy_boundaries,
                 const xt::xtensor<double, 1>& sig_t_r,
                 const xt::xtensor<double, 1>& sig_s_r, double awr_r);

  void add_background_nuclide(double sig_s, double awr);

  void solve();

  const xt::xtensor<double, 1>& energy_boundaries() const {
    return energy_boundaries_;
  }

  const xt::xtensor<double, 1>& avg_energy() const { return avg_energy_; }

  const xt::xtensor<double, 1>& sig_t() const { return sig_t_r_; }

  const xt::xtensor<double, 1>& sig_s() const { return sig_s_r_; }

  const xt::xtensor<double, 1>& flux() const { return flux_; }

  const xt::xtensor<double, 1>& chi() const { return chi_; }

  double awr() const { return A_r_; }

  double alpha() const { return a_r_; }

 private:
  xt::xtensor<double, 1> energy_boundaries_;
  xt::xtensor<double, 1> avg_energy_;
  xt::xtensor<double, 1> dlt_energy_;
  xt::xtensor<double, 1> flux_;
  xt::xtensor<double, 1> sig_t_r_;
  xt::xtensor<double, 1> sig_s_r_;
  xt::xtensor<double, 1> chi_;
  std::vector<BackgroundNuclide> background_nuclides_;
  double A_r_, a_r_;
  std::size_t g_max_r;

  double a_from_awr(double A) const {
    double res = (A - 1.) / (A + 1.);
    return res * res;
  }

  void generate_watt_spectrum();

  double compute_scatter_source(double Sgm1, std::size_t g);
};

}  // namespace scarabee

#endif
