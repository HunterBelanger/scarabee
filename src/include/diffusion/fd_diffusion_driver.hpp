#ifndef SCARABEE_FD_DIFFUSION_DRIVER_H
#define SCARABEE_FD_DIFFUSION_DRIVER_H

#include <diffusion_cross_section.hpp>
#include <diffusion/diffusion_geometry.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>

#include <memory>

namespace scarabee {

class FDDiffusionDriver {
  public:

  private:
    std::shared_ptr<DiffusionGeometry> geom_;
    xt::xtensor<double, 1> flux_; // Flux in each MAT tile in each group
    double keff_;
};

}

#endif