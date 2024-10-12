#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/pin_cell_type.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_PinCellType(py::module& m) {
  py::enum_<PinCellType>(m, "PinCellType")
      .value("Full", PinCellType::Full, "Full pin cell.")
      .value("XN", PinCellType::XN, "Half of a pin cell where x < 0.")
      .value("XP", PinCellType::XP, "Half of a pin cell where x > 0.")
      .value("YN", PinCellType::YN, "Half of a pin cell where y < 0.")
      .value("YP", PinCellType::YP, "Half of a pin cell where y > 0.")
      .value("I", PinCellType::I, "Quarter of a pin cell in quadrant I.")
      .value("II", PinCellType::II, "Quarter of a pin cell in quadrant II.")
      .value("III", PinCellType::III, "Quarter of a pin cell in quadrant III.")
      .value("IV", PinCellType::IV, "Quarter of a pin cell in quadrant IV.");
}