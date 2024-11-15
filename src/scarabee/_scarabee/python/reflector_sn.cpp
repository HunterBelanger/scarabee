#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <reflector_sn.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_ReflectorSN(py::module& m) {
  py::class_<ReflectorSN>(m, "ReflectorSN")
      .def(py::init<const std::vector<std::shared_ptr<CrossSection>>& /*xs*/,
                    const xt::xtensor<double, 1>& /*dx*/>())
      .def("solve", &ReflectorSN::solve)

      .def_property_readonly(
          "solved", &ReflectorSN::solved,
          "True if problem has been solved, False otherwsie.")

      .def_property_readonly("keff", &ReflectorSN::keff,
                             "Value of keff estimated by solver (1 by default "
                             "if no solution has been obtained).")

      .def_property_readonly("ngroups", &ReflectorSN::ngroups,
                             "Number of energy groups.")

      .def_property_readonly("size", &ReflectorSN::size, "Number of regions.")

      .def_property_readonly("nregions", &ReflectorSN::nregions,
                             "Number of regions.")

      .def_property_readonly("nsurfaces", &ReflectorSN::nsurfaces,
                             "Number of surfaces.")

      .def_property(
          "keff_tolerance", &ReflectorSN::keff_tolerance,
          &ReflectorSN::set_keff_tolerance,
          "Maximum relative absolute difference in keff for convergence")

      .def_property(
          "flux_tolerance", &ReflectorSN::flux_tolerance,
          &ReflectorSN::set_flux_tolerance,
          "Maximum relative absolute difference in flux for convergence")

      .def("flux", &ReflectorSN::flux,
           "Returns the scalar flux in group g in mesh region i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Flux in region i and in group g.\n",
           py::arg("i"), py::arg("g"))

      .def("current", &ReflectorSN::current,
           "Returns the net current in group g at surface i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Surface index.\n"
           "g : int\n"
           "    Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Net current at surface i and in group g.\n",
           py::arg("i"), py::arg("g"))

      .def("volume", &ReflectorSN::volume,
           "Returns the volume of region i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Volume of region i.\n",
           py::arg("i"))

      .def("xs", &ReflectorSN::xs,
           "Returns the CrossSection in region i.\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "    Region index.\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "    Material cross sections at r.\n",
           py::arg("i"))

      .def("homogenize", &ReflectorSN::homogenize,
           "Computes a homogenized set of cross sections for the set of "
           "provided region indices.\n\n"
           "Parameters\n"
           "----------\n"
           "regions : list of int\n"
           "          List of regions for homogenization.\n"
           "Returns\n"
           "-------\n"
           ":py:class:`CrossSection`\n"
           "        Homogenized cross section.\n",
           py::arg("regions"))

      .def("homogenize_flux_spectrum", &ReflectorSN::homogenize_flux_spectrum,
           "Computes a homogenized flux spectrum based on the list of "
           "provided region indices. This method will raise an exception if "
           "the problem has not yet been solved.\n\n"
           "Parameters\n"
           "----------\n"
           "regions : list of int\n"
           "          List of regions for homogenization.\n"
           "Returns\n"
           "-------\n"
           "ndarray of floats\n"
           "                 Homogenized flux spectrum.",
           py::arg("regions"));
}