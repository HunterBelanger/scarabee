#ifndef SCARABEE_CONSTANTS_H
#define SCARABEE_CONSTANTS_H

#include <limits>

namespace scarabee {

constexpr double PI{3.14159265358979323846264338327950288};
constexpr double PI_2{0.5 * PI};
constexpr double INF{std::numeric_limits<double>::max()};
constexpr double SURFACE_COINCIDENT{1E-12};
constexpr double VEC_FP_TOL{1.E-10};
constexpr std::size_t MAX_SURFS{5};

}  // namespace scarabee

#endif
