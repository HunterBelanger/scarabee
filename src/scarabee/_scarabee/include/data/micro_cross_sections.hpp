#ifndef SCARABEE_MICRO_CROSS_SECTIONS_H
#define SCARABEE_MICRO_CROSS_SECTIONS_H

#include <data/xs1d.hpp>

#include <xtensor/xtensor.hpp>

#include <optional>

namespace scarabee {

struct MicroDepletionXS {
  std::optional<XS1D> n_fission;
  std::optional<XS1D> n_gamma;
  std::optional<XS1D> n_2n;
  std::optional<XS1D> n_3n;
  std::optional<XS1D> n_a;
  std::optional<XS1D> n_p;
};

struct MicroNuclideXS {
  XS1D Et;
  XS1D Dtr;
  XS2D Es;
  XS1D Ea;
  XS1D Ef;
  XS1D nu;
  XS1D chi;
};

struct ResonantOneGroupXS {
  double Dtr;
  double Ea;
  double Ef;
  xt::xtensor<double, 2> Es; // First index is legendre moment, second is outgoing energy group
  std::size_t gout_min; // Index of first tabulated outgoing energy group

  std::optional<double> n_gamma;

  // Scarabée assumes that (n,2n), (n,3n), (n,a), and (n,p) are not resonant
  // i.e. not dilution dependent.
};

}

#endif