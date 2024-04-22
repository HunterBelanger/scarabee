#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/simple_pin_cell.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_SimplePinCell(py::module& m) {
  py::class_<SimplePinCell, Cell, std::shared_ptr<SimplePinCell>>(m, "SimplePinCell")
      .def(py::init<const std::vector<double>& /*rads*/,
                    const std::vector<std::shared_ptr<TransportXS>>& /*mats*/,
                    double /*dx*/, double /*dy*/>(),
           "An annular pin cell centered in the cell with no angular segments.\n"
           "Must provide one more TransportXS than radii, as that material\n"
           "will fill the cell out to the boundary.\n\n"
           "Arguments:\n"
           "    radii  List of radii of annular regions\n"
           "    mats   List of TransportXS\n"
           "    dx     Width of cell in x\n"
           "    dy     Height of cell in y",
           py::arg("radii"), py::arg("mats"), py::arg("dx"), py::arg("dy"));
}