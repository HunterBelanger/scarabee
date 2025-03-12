#ifndef PIN_CELL_TYPE_H
#define PIN_CELL_TYPE_H

#include <cstdint>

namespace scarabee {

enum class PinCellType : std::uint8_t { Full, XP, XN, YP, YN, I, II, III, IV };

enum class BWRCornerType : std::uint8_t { I, II, III, IV };

}  // namespace scarabee

#endif