#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <diffusion/fd_diffusion_driver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_FDDiffusionDriver(py::module& m) {
  py::class_<FDDiffusionDriver>(m, "FDDiffusionDriver",
  "A FDDiffusionDriver solved a diffusion problem using the cell centered "
  "finite-difference method. It is capable of solving 1D, 2D, or 3D problems "
  "which are dfined by providing a :py:class:`DiffusionGeometry` instance.")
  
  .def(py::init<std::shared_ptr<DiffusionGeometry> /*geom*/>(),
  "Initializes a finite-difference diffusion solver.\n\n"
  "Parameters\n"
  "----------\n"
  "geom : DiffusionGeometry\n"
  "       Problem deffinition to solve.")
  
  .def("solve", &FDDiffusionDriver::solve, "Solves the diffusion problem.")
  
  .def_property_readonly("ngroups", &FDDiffusionDriver::ngroups, "Number of energy groups.")
  
  .def_property_readonly("solved", &FDDiffusionDriver::solved, "True if the problem has been solved, False otherwise.")
  
  .def_property_readonly("keff", &FDDiffusionDriver::keff, "Value of keff. This is 1 by default is solved is False.")

  .def_property_readonly("flux", &FDDiffusionDriver::flux)
  
  .def_property("keff_tolerance", &FDDiffusionDriver::keff_tolerance, &FDDiffusionDriver::set_keff_tolerance, "Maximum relative error in keff for problem convergence.")
  
  .def_property("flux_tolerance", &FDDiffusionDriver::flux_tolerance, &FDDiffusionDriver::flux_tolerance, "Maximum relative error in the flux for problem convergence.");
}