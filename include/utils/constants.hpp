#ifndef SCARABEE_CONSTANTS_H
#define SCARABEE_CONSTANTS_H

#include <cstdint>
#include <limits>

constexpr double PI{3.14159265358979323846264338327950288};
constexpr double PI_2{0.5 * PI};
constexpr double INF{std::numeric_limits<double>::max()};
constexpr double SURFACE_COINCIDENT{1E-12};
constexpr std::size_t MAX_SURFS{5};

#endif
