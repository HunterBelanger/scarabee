#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <cylindrical_flux_solver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CylindricalFluxSolver(py::module& m) {
  py::class_<CylindricalFluxSolver, std::shared_ptr<CylindricalFluxSolver>>(
      m, "CylindricalFluxSolver")
      .def(py::init<std::shared_ptr<CylindricalCell>>(),
           "Flux solver for a 1D annular cell problem.\n\n"
           "Arguments:\n"
           "    cell  CylindricalCell for which to find the flux",
           py::arg("cell"))

      .def_property_readonly("ngroups", &CylindricalFluxSolver::ngroups,
           "Number of energy groups")

      .def_property_readonly("nregions", &CylindricalFluxSolver::nregions,
           "Number of annular regions")

      .def("flux", &CylindricalFluxSolver::flux,
           "Flux in group g region i.\n\n"
           "Arguments:\n"
           "    g  group index\n"
           "    i  region index",
           py::arg("g"), py::arg("i"))

      .def_property("flux_tolerance", &CylindricalFluxSolver::flux_tolerance, &CylindricalFluxSolver::set_flux_tolerance,
           "Tolerance for flux convergence")

      .def("set_flux_tolerance", &CylindricalFluxSolver::set_flux_tolerance,
           "Sets the flux convervence tolerance.\n\n"
           "Arguments:\n"
           "    flux_tol  desired absolute relative error in flux",
           py::arg("flux_tol"))

      .def_property_readonly("keff", &CylindricalFluxSolver::keff,
           "Multiplication factor of cell")

      .def_property("keff_tolerance", &CylindricalFluxSolver::keff_tolerance, &CylindricalFluxSolver::set_keff_tolerance,
           "Tolerance for keff convergence")

      .def_property("albedo", &CylindricalFluxSolver::albedo, &CylindricalFluxSolver::set_albedo, "Albedo for outer cell boundary")

      .def("j_ext", &CylindricalFluxSolver::j_ext,
           "External current.\n\n"
           "Arguments:\n"
           "    g  group index",
           py::arg("g"))

      .def("set_j_ext", &CylindricalFluxSolver::set_j_ext,
           "Sets the external current.\n\n"
           "Arguments:\n"
           "    g  group index\n"
           "    j  external current",
           py::arg("g"), py::arg("j"))

      .def("j_neg", &CylindricalFluxSolver::j_neg,
           "Inward current into cell.\n\n"
           "Arguments:\n"
           "    g  group index",
           py::arg("g"))

      .def("j_pos", &CylindricalFluxSolver::j_pos,
           "Outward current from cell.\n\n"
           "Arguments:\n"
           "    g  group index",
           py::arg("g"))

      .def("j", &CylindricalFluxSolver::j,
           "Net current out of cell.\n\n"
           "Arguments:\n"
           "    g  group index",
           py::arg("g"))

      .def("solve", &CylindricalFluxSolver::solve,
           "Solves the system for the flux and eigenvalue.")

      .def_property_readonly("solved", &CylindricalFluxSolver::solved,
           "Returns true if the system has been solved.");
}
