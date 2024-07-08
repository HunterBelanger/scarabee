#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <diffusion/fd_diffusion_driver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_FDDiffusionDriver(py::module& m) {
  py::class_<FDDiffusionDriver>(
      m, "FDDiffusionDriver",
      "A FDDiffusionDriver solves a diffusion problem using the cell centered "
      "finite-difference method. It is capable of solving 1D, 2D, or 3D "
      "problems which are defined by providing a "
      ":py:class:`DiffusionGeometry` instance.")

      .def(py::init<std::shared_ptr<DiffusionGeometry> /*geom*/>(),
           "Initializes a finite-difference diffusion solver.\n\n"
           "Parameters\n"
           "----------\n"
           "geom : DiffusionGeometry\n"
           "       Problem deffinition to solve.")

      .def("solve", &FDDiffusionDriver::solve, "Solves the diffusion problem.")

      .def_property_readonly("ngroups", &FDDiffusionDriver::ngroups,
                             "Number of energy groups.")

      .def_property_readonly(
          "solved", &FDDiffusionDriver::solved,
          "True if the problem has been solved, False otherwise.")

      .def_property_readonly(
          "keff", &FDDiffusionDriver::keff,
          "Value of keff. This is 1 by default is solved is False.")

      .def_property("keff_tolerance", &FDDiffusionDriver::keff_tolerance,
                    &FDDiffusionDriver::set_keff_tolerance,
                    "Maximum relative error in keff for problem convergence.")

      .def_property(
          "flux_tolerance", &FDDiffusionDriver::flux_tolerance,
          &FDDiffusionDriver::flux_tolerance,
          "Maximum relative error in the flux for problem convergence.")

      .def("flux", &FDDiffusionDriver::flux,
           "Returns the computed flux, along with the mesh bounds. The first "
           "dimension "
           "of the flux array is the energy group index. The second index is "
           "for the x "
           "coordinate. If the problem is 2 or 3 dimensional, the third and "
           "fourth "
           "indices are for the y and z coordinates respectively.\n\n"
           "Returns\n"
           "-------\n"
           "flux : ndarray\n"
           "       2D, 3D, or 4D array containing the multigroup flux.\n"
           "x_bounds : ndarray\n"
           "           1D array with the x-bounds for the flux mesh.\n"
           "y_bounds : ndarray or None\n"
           "           1D array with the y-bounds for the flux mesh, if a 2D "
           "problem.\n"
           "z_bounds : ndarray or None\n"
           "           1D array with the z-bounds for the flux mesh, if a 3D "
           "problem.\n")

      .def("power", &FDDiffusionDriver::power,
           "Returns the computed power distribution, along with the mesh "
           "bounds. The first dimension of the power array is for the x "
           "coordinate. If the problem is 2 or 3 dimensional, the second and "
           "third indices are for the y and z coordinates respectively.\n\n"
           "Returns\n"
           "-------\n"
           "power : ndarray\n"
           "        1D, 2D, or 3D array containing the power distribution.\n"
           "x_bounds : ndarray\n"
           "           1D array with the x-bounds for the power mesh.\n"
           "y_bounds : ndarray or None\n"
           "           1D array with the y-bounds for the power mesh, if a 2D "
           "problem.\n"
           "z_bounds : ndarray or None\n"
           "           1D array with the z-bounds for the power mesh, if a 3D "
           "problem.\n");
}
