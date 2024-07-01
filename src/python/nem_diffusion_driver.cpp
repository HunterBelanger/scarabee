#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <diffusion/nem_diffusion_driver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_NEMDiffusionDriver(py::module& m) {
  py::class_<NEMDiffusionDriver>(
      m, "NEMDiffusionDriver",
      "A NEMDiffusionDriver solves a diffusion problem using the nodal "
      "expansion method. It is capable of solving 3D problems which are "
      "defined by providing a :py:class:`DiffusionGeometry` instance.")

      .def(py::init<std::shared_ptr<DiffusionGeometry> /*geom*/>(),
           "Initializes a nodal diffusion solver.\n\n"
           "Parameters\n"
           "----------\n"
           "geom : DiffusionGeometry\n"
           "       Problem deffinition to solve.")

      .def("solve", &NEMDiffusionDriver::solve, "Solves the diffusion problem.")

      .def_property_readonly("ngroups", &NEMDiffusionDriver::ngroups,
                             "Number of energy groups.")

      .def_property_readonly(
          "solved", &NEMDiffusionDriver::solved,
          "True if the problem has been solved, False otherwise.")

      .def_property_readonly(
          "keff", &NEMDiffusionDriver::keff,
          "Value of keff. This is 1 by default is solved is False.")

      .def_property("keff_tolerance", &NEMDiffusionDriver::keff_tolerance,
                    &NEMDiffusionDriver::set_keff_tolerance,
                    "Maximum relative error in keff for problem convergence.")

      .def_property(
          "flux_tolerance", &NEMDiffusionDriver::flux_tolerance,
          &NEMDiffusionDriver::flux_tolerance,
          "Maximum relative error in the flux for problem convergence.");
}
