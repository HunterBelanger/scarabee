#ifndef REFLECTOR_SN_H
#define REFLECTOR_SN_H

#include <cross_section.hpp>

#include <xtensor/xtensor.hpp>

#include <array>
#include <memory>

namespace scarabee {

class ReflectorSN {
  public:
    ReflectorSN(const std::vector<std::shared_ptr<CrossSection>>& xs,
                const xt::xtensor<double, 1>& dx):
                xs_(xs), dx_(dx) {}

    void solve();

  private:
    std::vector<std::shared_ptr<CrossSection>> xs_;
    xt::xtensor<double, 1> dx_;
    xt::xtensor<double, 2> flux_; // group, spatial bin
    xt::xtensor<double, 2> Q_;    // group, spatial bin
    double keff_;


    std::array<double, 16> mu_ {-0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374, 0.0950125098376374,
    0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 0.7554044083550030, 0.8656312023878318,
    0.9445750230732326, 0.9894009349916499};

    std::array<double, 16> wgt_ {0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
    0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685, 0.1894506104550685,
    0.1826034150449236, 0.1691565193950025, 0.1495959888165767, 0.1246289712555339, 0.0951585116824928,
    0.0622535239386479, 0.0271524594117541};


    void sweep(xt::xtensor<double, 2>& flux, const xt::xtensor<double, 2>& Q);
    double calc_keff(const xt::xtensor<double, 2>& old_flux, const xt::xtensor<double, 2>& new_flux, const double keff) const;
    void fill_fission_source(xt::xtensor<double, 2>& Qfiss, const xt::xtensor<double, 2>& flux) const;
    void fill_scatter_source(xt::xtensor<double, 2>& Qscat, const xt::xtensor<double, 2>& flux) const;
};

}

#endif