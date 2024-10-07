#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/pin_cell_type.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_PinCellType(py::module& m) {
  py::enum_<PinCellType>(m, "PinCellType")
      .value("Full", PinCellType::Full)
      .value("XN", PinCellType::XN)
      .value("XP", PinCellType::XP)
      .value("YN", PinCellType::YN)
      .value("YP", PinCellType::YP)
      .value("I", PinCellType::I)
      .value("II", PinCellType::II)
      .value("III", PinCellType::III)
      .value("IV", PinCellType::IV);
}