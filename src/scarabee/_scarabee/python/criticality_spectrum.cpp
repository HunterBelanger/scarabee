#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include <utils/criticality_spectrum.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CriticalitySpectrum(py::module& m) {
  py::class_<CriticalitySpectrum>(m, "CriticalitySpectrum")
      .def_property_readonly("ngroups", &CriticalitySpectrum::ngroups,
                             "Number of energy groups")

      .def_property_readonly("k_inf", &CriticalitySpectrum::k_inf,
                             "Infinite multiplication factor")

      .def_property_readonly("B2", &CriticalitySpectrum::B2,
                             "Critical buckling B^2")

      .def_property_readonly("buckling", &CriticalitySpectrum::buckling,
                             "Critical buckling B^2")

      .def_property_readonly(
          "flux", py::overload_cast<>(&CriticalitySpectrum::flux, py::const_),
          py::return_value_policy::reference_internal,
          "Array contianing the flux spectrum")

      .def_property_readonly(
          "current",
          py::overload_cast<>(&CriticalitySpectrum::current, py::const_),
          py::return_value_policy::reference_internal,
          "Array contianing the current spectrum")

      .def_property_readonly(
          "diff_coeff",
          py::overload_cast<>(&CriticalitySpectrum::diff_coeff, py::const_),
          py::return_value_policy::reference_internal,
          "Array contianing the diffusion coefficients");

  py::class_<P1CriticalitySpectrum, CriticalitySpectrum>(
      m, "P1CriticalitySpectrum")
      .def(py::init<std::shared_ptr<CrossSection>>(),
           "Computes the criticality energy spectrum using the P1 leakage "
           "approximation.\n\n"
           "Parameters\n"
           "----------\n"
           "xs : CrossSection\n"
           "     Homogenized set of cross sections for the system.\n",
           py::arg("xs"));

  py::class_<B1CriticalitySpectrum, CriticalitySpectrum>(
      m, "B1CriticalitySpectrum")
      .def(py::init<std::shared_ptr<CrossSection>>(),
           "Computes the criticality energy spectrum using the B1 leakage "
           "approximation.\n\n"
           "Parameters\n"
           "----------\n"
           "xs : CrossSection\n"
           "     Homogenized set of cross sections for the system.\n",
           py::arg("xs"));
}
