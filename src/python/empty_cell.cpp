#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/empty_cell.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_EmptyCell(py::module& m) {
  py::class_<EmptyCell, Cell, std::shared_ptr<EmptyCell>>(m, "EmptyCell")
      .def(py::init<const std::shared_ptr<CrossSection>& /*mat*/, double /*dx*/,
                    double /*dy*/>(),
           "An empty cell with one flat source region.\n\n"
           "Arguments:\n"
           "    mat    CrossSection for cell\n"
           "    dx     Width of cell in x\n"
           "    dy     Height of cell in y",
           py::arg("mat"), py::arg("dx"), py::arg("dy"));
}
