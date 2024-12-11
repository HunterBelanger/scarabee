#ifndef SCARABEE_WATER_H
#define SCARABEE_WATER_H

#include <data/material.hpp>

#include <memory>
#include <optional>

namespace scarabee {

double water_density(double temperature, double pressure);

std::shared_ptr<Material> borated_water(double boron_ppm, double temperature,
                                        double pressure,
                                        std::shared_ptr<NDLibrary> ndl);

}  // namespace scarabee

#endif