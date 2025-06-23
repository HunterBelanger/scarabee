#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/cell.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_Cell(py::module& m) {
  py::class_<Cell, std::shared_ptr<Cell>>(m, "Cell")

      .def("inside", &Cell::inside,
           "Checks if a position - direction pair are inside the cell.\n\n"
           "Parameters\n"
           "----------\n"
           "r : Vector\n"
           "    Position to test.\n"
           "u : Direction\n"
           "    Direction vector for disambiguating a the region.\n\n"
           "Returns\n"
           "-------\n"
           "bool\n"
           "     True is r and u are in the cell, False otherwise.\n",
           py::arg("r"), py::arg("u"))

      .def("distance", &Cell::distance,
           "Distance that can be traveled within cell for given position and "
           "direction.\n\n"
           "Parameters\n"
           "----------\n"
           "r : Vector\n"
           "    Starting position.\n"
           "u : Direction\n"
           "    Direction of travel in the cell.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Distance that can be traveled.\n",
           py::arg("r"), py::arg("u"))

      .def("get_all_fsr_ids", &Cell::get_all_fsr_ids,
           "Returns a set containing the IDs of all flat source regions in the "
           "cell.\n\n"
           "Returns\n"
           "-------\n"
           "set of int\n"
           "    IDs of all FSRs in the cell.")

      .def_property_readonly("num_fsrs", &Cell::num_fsrs,
                             "Number of flat source regions in the cell.")

      .def_property_readonly("dx", &Cell::dx, "Width of cell along x.")

      .def_property_readonly("dy", &Cell::dy, "Width of cell along y.");
}
