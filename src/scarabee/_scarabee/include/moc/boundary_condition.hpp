#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <cstdint>

namespace scarabee {

enum class BoundaryCondition : std::uint8_t { Reflective, Vacuum, Periodic };

}

#endif
