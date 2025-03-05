#ifndef SCARABEE_MICRO_CROSS_SECTIONS_H
#define SCARABEE_MICRO_CROSS_SECTIONS_H

#include <data/xs1d.hpp>
#include <data/xs2d.hpp>
#include <utils/serialization.hpp>

#include <xtensor/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/optional.hpp>

#include <optional>

namespace scarabee {

struct MicroDepletionXS {
  std::optional<XS1D> n_fission{std::nullopt};
  std::optional<XS1D> n_gamma{std::nullopt};
  std::optional<XS1D> n_2n{std::nullopt};
  std::optional<XS1D> n_3n{std::nullopt};
  std::optional<XS1D> n_a{std::nullopt};
  std::optional<XS1D> n_p{std::nullopt};

  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(n_fission), CEREAL_NVP(n_gamma), CEREAL_NVP(n_2n),
        CEREAL_NVP(n_3n), CEREAL_NVP(n_a), CEREAL_NVP(n_p));
  }
};

struct MicroNuclideXS {
  XS1D Et;
  XS1D Dtr;
  XS2D Es;
  XS1D Ea;
  XS1D Ef;
  XS1D nu;
  XS1D chi;

  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(Et), CEREAL_NVP(Dtr), CEREAL_NVP(Es), CEREAL_NVP(Es),
        CEREAL_NVP(Ef), CEREAL_NVP(nu), CEREAL_NVP(chi));
  }
};

struct ResonantOneGroupXS {
  double Dtr{0.};
  double Ea{0.};
  double Ef{0.};
  xt::xtensor<double, 2>
      Es;  // First index is legendre moment, second is outgoing energy group
  std::size_t gout_min{0};  // Index of first tabulated outgoing energy group

  std::optional<double> n_gamma{std::nullopt};

  // Scarabée assumes that (n,2n), (n,3n), (n,a), and (n,p) are not resonant
  // i.e. not dilution dependent.
};

}  // namespace scarabee

#endif