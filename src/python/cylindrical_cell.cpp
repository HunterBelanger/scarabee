#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <cylindrical_cell.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CylindricalCell(py::module& m) {
  py::class_<CylindricalCell, std::shared_ptr<CylindricalCell>>(m, "CylindricalCell")
      .def(
          py::init<const std::vector<double>& /*radii*/,
                   const std::vector<std::shared_ptr<TransportXS>>& /*mats*/>(),
          "1D annular cell problem.\n\n"
          "Arguments:\n"
          "    radii  List of region radii\n"
          "    mats   List of TransportXS for each region.",
          py::arg("radii"), py::arg("mats"))

      .def("solved", &CylindricalCell::solved,
           "Returns true is the system has been solved")

      .def("solve", &CylindricalCell::solve,
           "Solves the system for partial flux responses")

      .def("ngroups", &CylindricalCell::ngroups, "Number of energy groups")

      .def("nregions", &CylindricalCell::nregions, "Number of annular regions")

      .def("Sb", &CylindricalCell::Sb,
           "Outer suurface area (circumfrence) of cell")

      .def("V", &CylindricalCell::V,
           "Volume of region i\n\n"
           "Arguments:\n"
           "    i  region index",
           py::arg("i"))

      .def("R", &CylindricalCell::R,
           "Radius of region i\n\n"
           "Arguments:\n"
           "    i  region index",
           py::arg("i"))

      .def("Y", &CylindricalCell::Y,
           "Current response flux. This is the flux in region i\n"
           "cause by a neutron entering the cell from the outer\n"
           "surface (with albedo a) in group g.\n\n"
           "Arguments:\n"
           "    a  albedo of boundary\n"
           "    g  group index\n"
           "    i  region index",
           py::arg("a"), py::arg("g"), py::arg("i"))

      .def("x", &CylindricalCell::x,
           "First escape to the boundary from region i.\n\n"
           "Arguments:\n"
           "    g  group index\n"
           "    i  region index",
           py::arg("g"), py::arg("i"))

      .def("X", &CylindricalCell::X,
           "Source response flux. This is the flux in region i\n"
           "cause by a unit source in region k, in group g.\n\n"
           "Arguments:\n"
           "    a  albedo of boundary\n"
           "    g  group index\n"
           "    i  region index for response\n"
           "    k  region index for source",
           py::arg("a"), py::arg("g"), py::arg("i"), py::arg("k"))

      .def("Gamma", &CylindricalCell::Gamma,
           "Multicollision blackness. The probability of neutrons\n"
           "entering through the outer surface Sb are absorbed.\n\n"
           "Arguments:\n"
           "    g  group index",
           py::arg("g"))

      .def("p", &CylindricalCell::p,
           "Symmetric collision probability.\n\n"
           "Arguments:\n"
           "    g  group index\n"
           "    i  origin region index\n"
           "    k  final region index",
           py::arg("g"), py::arg("i"), py::arg("k"))

      .def("mat", &CylindricalCell::mat,
           "TransportXS for region i.\n\n"
           "Arguments:\n"
           "    i  region index",
           py::arg("i"));
}
