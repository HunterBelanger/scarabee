#include <pybind11/pybind11.h>

#include <utils/version.hpp>

namespace py = pybind11;

extern void init_TransportXS(py::module&);
extern void init_CylindricalCell(py::module&);
extern void init_CylindricalFluxSolver(py::module&);
extern void init_Logging(py::module&);
extern void init_Vector(py::module&);

PYBIND11_MODULE(pyScarabee, m) {
  init_Logging(m);
  init_Vector(m);
  init_TransportXS(m);
  init_CylindricalCell(m);
  init_CylindricalFluxSolver(m);

  m.attr("__author__") = "Hunter Belanger";
  m.attr("__copyright__") = "Copyright 2024, Hunter Belanger";
  m.attr("__license__") = "GPL-3.0-or-later";
  m.attr("__maintainer__") = "Hunter Belanger";
  m.attr("__email__") = "hunter.belanger@gmail.com";
  m.attr("__version__") = scarabee::VERSION_STRING;
}
