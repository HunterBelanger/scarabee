#ifndef SCARABEE_CRITICALITY_SPECTRUM_H
#define SCARABEE_CRITICALITY_SPECTRUM_H

#include <cross_section.hpp>

#include <xtensor/xtensor.hpp>

namespace scarabee {

xt::xtensor<double, 2> P1_spectrum(std::shared_ptr<CrossSection> xs);

}

#endif