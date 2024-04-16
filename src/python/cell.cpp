#include <pybind11/pybind11.h>

#include <moc/cell.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_Cell(py::module& m) {
  py::class_<Cell, std::shared_ptr<Cell>>(m, "Cell")

      .def("inside", &Cell::inside,
           "Checks if a position - direction pair are inside the cell.\n\n"
           "Arguments:\n"
           "    r  Vector of position\n"
           "    u  Direction at r for boundary",
           py::arg("r"), py::arg("u"))

      .def("distance", &Cell::distance,
           "Distance within cell for given position and direction.\n\n"
           "Arguments:\n"
           "    r  starting position\n"
           "    u  direction of travel",
           py::arg("r"), py::arg("u"))

      .def("dx", &Cell::dx, "Width of Cell in x.")

      .def("dy", &Cell::dy, "Width of Cell in y.");
}
