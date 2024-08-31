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
                const xt::xtensor<double, 1>& dx);

    void solve();
    bool solved() const { return solved_; }

    double keff() const { return keff_; }

    double keff_tolerance() const { return keff_tol_; }
    void set_keff_tolerance(double ktol);

    double flux_tolerance() const { return flux_tol_; }
    void set_flux_tolerance(double ftol);

    std::size_t size() const { return xs_.size(); }
    std::size_t nregions() const { return xs_.size(); }
    std::size_t ngroups() const { return ngroups_; }

    const std::shared_ptr<CrossSection> xs(std::size_t i) const; 
    double volume(std::size_t i) const;
    
    double flux(std::size_t i, std::size_t g) const;

    std::shared_ptr<CrossSection> homogenize(const std::vector<std::size_t>& regions) const;
    xt::xtensor<double, 1> homogenize_flux_spectrum(const std::vector<std::size_t>& regions) const;

  private:
    std::vector<std::shared_ptr<CrossSection>> xs_;
    xt::xtensor<double, 1> dx_;
    xt::xtensor<double, 2> flux_; // group, spatial bin
    xt::xtensor<double, 2> Q_;    // group, spatial bin
    double keff_ {1.}; 
    double keff_tol_ {1.E-5};
    double flux_tol_ {1.E-5};
    std::size_t ngroups_;
    bool solved_ {false};

    void sweep(xt::xtensor<double, 2>& flux, const xt::xtensor<double, 2>& Q);
    double calc_keff(const xt::xtensor<double, 2>& old_flux, const xt::xtensor<double, 2>& new_flux, const double keff) const;
    void fill_fission_source(xt::xtensor<double, 2>& Qfiss, const xt::xtensor<double, 2>& flux) const;
    void fill_scatter_source(xt::xtensor<double, 2>& Qscat, const xt::xtensor<double, 2>& flux) const;

    static const std::array<double, 32> mu_;
    static const std::array<double, 32> wgt_;
};

}

#endif