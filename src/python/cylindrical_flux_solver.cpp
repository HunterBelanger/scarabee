#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <cylindrical_flux_solver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CylindricalFluxSolver(py::module& m) {
  py::class_<CylindricalFluxSolver, std::shared_ptr<CylindricalFluxSolver>>(
      m, "CylindricalFluxSolver",
      "A CylindricalFluxSolver object computes the static flux in a 1D annular "
      "pin cell, based on the collision probabilities and response fluxes "
      "provided by a :py:class:`CylindricalCell` instance.")

      .def(py::init<std::shared_ptr<CylindricalCell>>(),
           "Flux solver for a 1D annular cell problem.\n\n"
           "Parameters\n"
           "----------\n"
           "cell : CylindricalCell\n"
           "       Pin cell in which to find the flux.",
           py::arg("cell"))

      .def_property_readonly("ngroups", &CylindricalFluxSolver::ngroups,
                             "Number of energy groups.")

      .def_property_readonly("nregions", &CylindricalFluxSolver::nregions,
                             "Number of annular regions.")

      .def("flux", &CylindricalFluxSolver::flux,
           "Flux in region i, group g.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Flux in region i group g.",
           py::arg("i"), py::arg("g"))

      .def("volume", &CylindricalFluxSolver::volume,
           "Volume of region i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Volume of region i.",
           py::arg("i"))

      .def("xs", &CylindricalFluxSolver::xs,
           "CrossSection in region i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "            Cross sections in region i.",
           py::arg("i"))

      .def_property("flux_tolerance", &CylindricalFluxSolver::flux_tolerance,
                    &CylindricalFluxSolver::set_flux_tolerance,
                    "Tolerance for flux convergence.")

      .def_property_readonly(
          "keff", &CylindricalFluxSolver::keff,
          "Multiplication factor of cell (only valid if solved is True).")

      .def_property("keff_tolerance", &CylindricalFluxSolver::keff_tolerance,
                    &CylindricalFluxSolver::set_keff_tolerance,
                    "Tolerance for keff convergence.")

      .def_property("albedo", &CylindricalFluxSolver::albedo,
                    &CylindricalFluxSolver::set_albedo,
                    "Albedo for outer cell boundary.")

      .def("j_ext", &CylindricalFluxSolver::j_ext,
           "External current at boundary in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     External current in group g.",
           py::arg("g"))

      .def("set_j_ext", &CylindricalFluxSolver::set_j_ext,
           "Sets the external current at the boundary in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group index.\n"
           "j : float\n"
           "    External current.",
           py::arg("g"), py::arg("j"))

      .def("j_neg", &CylindricalFluxSolver::j_neg,
           "Inward current into the cell in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Inward current in group g.",
           py::arg("g"))

      .def("j_pos", &CylindricalFluxSolver::j_pos,
           "Outward current from the cell in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Outward current in group g.",
           py::arg("g"))

      .def("j", &CylindricalFluxSolver::j,
           "Net current out of cell in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Net current in group g.",
           py::arg("g"))

      .def("solve", &CylindricalFluxSolver::solve,
           "Solves the system for the flux and eigenvalue.")

      .def("homogenize", &CylindricalFluxSolver::homogenize,
           "Computes a homogenized set of cross sections for the problem.\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "            Homogenized cross section.")

      .def("homogenize_flux_spectrum",
           &CylindricalFluxSolver::homogenize_flux_spectrum,
           "Computes a homogenized flux spectrum which can be used for energy "
           "condensation.\n\n"
           "Returns\n"
           "-------\n"
           "list of floats\n"
           "               Homogenized flux spectrum.")

      .def_property_readonly(
          "solved", &CylindricalFluxSolver::solved,
          "True if the system has been solved, False otherwise.");
}
