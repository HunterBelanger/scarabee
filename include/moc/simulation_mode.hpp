#ifndef SCARABEE_SIMULATION_MODE_H
#define SCARABEE_SIMULATION_MODE_H

#include <cstdint>

namespace scarabee
{

enum class SimulationMode : std::uint8_t { FixedSource, Keff };
    
} // namespace scarabee


#endif