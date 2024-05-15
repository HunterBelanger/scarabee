#ifndef SCARABEE_ND_LIBRARY_H
#define SCARABEE_ND_LIBRARY_H

#include <data/material.hpp>
#include <cross_section.hpp>

#include <xtensor/xtensor.hpp>
#include <highfive/highfive.hpp>

namespace H5 = HighFive;

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace scarabee {

class NDLibrary;

struct NuclideHandle {
  std::string name;
  std::string label;
  std::vector<double> temperatures;
  std::vector<double> dilutions;
  double awr;
  double potential_xs;
  std::uint32_t ZA;
  bool fissile;
  bool resonant;

  std::shared_ptr<xt::xtensor<double, 3>> absorption;
  std::shared_ptr<xt::xtensor<double, 4>> scatter;
  std::shared_ptr<xt::xtensor<double, 4>> p1_scatter;
  std::shared_ptr<xt::xtensor<double, 3>> fission;
  std::shared_ptr<xt::xtensor<double, 2>> chi;
  std::shared_ptr<xt::xtensor<double, 2>> nu;

  bool loaded() const { return absorption != nullptr; }
  void load_xs_from_hdf5(const NDLibrary& ndl);
};

class NDLibrary {
 public:
  NDLibrary(const std::string& fname);

  std::size_t ngroups() const { return ngroups_; }

  const std::string& library() const { return library_; }

  const std::string& group_structure() const { return group_structure_; }

  const std::vector<double>& group_bounds() const { return group_bounds_; }

  NuclideHandle& get_nuclide(const std::string& name);
  const NuclideHandle& get_nuclide(const std::string& name) const;

  void convert_fractions(Material& mat) const;

  std::shared_ptr<CrossSection> interp_nuclide_xs(const std::string& name,
                                                  const double temp,
                                                  const double dil);

  std::shared_ptr<CrossSection> carlvik_two_term(const std::string& name,
                                                 const double mat_pot_xs,
                                                 const double temp,
                                                 const double N, const double C,
                                                 const double Ee);

  const std::shared_ptr<H5::File>& h5() const { return h5_; }

 private:
  std::map<std::string, NuclideHandle> nuclide_handles_;
  std::vector<double> group_bounds_;
  std::string library_;
  std::string group_structure_;
  std::size_t ngroups_;
  std::shared_ptr<H5::File> h5_;

  NDLibrary(const NDLibrary&) = delete;
  NDLibrary& operator=(const NDLibrary&) = delete;

  void get_temp_interp_params(double temp, const NuclideHandle& nuc,
                              std::size_t& i, double& f) const;
  void get_dil_interp_params(double dil, const NuclideHandle& nuc,
                             std::size_t& i, double& f) const;

  void interp_1d(xt::xtensor<double, 1>& E, const xt::xtensor<double, 2> nE,
                 std::size_t it, double f_temp) const;
  void interp_1d(xt::xtensor<double, 1>& E, const xt::xtensor<double, 3> nE,
                 std::size_t it, double f_temp, std::size_t id,
                 double f_dil) const;
  void interp_2d(xt::xtensor<double, 2>& E, const xt::xtensor<double, 4> nE,
                 std::size_t it, double f_temp, std::size_t id,
                 double f_dil) const;
};

}  // namespace scarabee

#endif