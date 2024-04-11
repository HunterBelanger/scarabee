#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include <transport_xs.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_TransportXS(py::module& m) {
  py::class_<TransportXS, std::shared_ptr<TransportXS>>(m, "TransportXS")
      .def(py::init<const xt::xtensor<double, 1>& /*Et*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const xt::xtensor<double, 1>& /*vEf*/,
                    const xt::xtensor<double, 1>& /*chi*/>(),
           "Transport corrected macroscopic cross sections.\n\n"
           "Arguments:\n"
           "    Et    transport corrected total cross section\n"
           "    Ea    absorption cross section\n"
           "    Es    transport corrected scattering cross section matrix\n"
           "    vEf   fission yield time cross section\n"
           "    chi   fission spectrum",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"), py::arg("vEf"),
           py::arg("chi"))

      .def(py::init<const xt::xtensor<double, 1>& /*Et*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/>(),
           "Transport corrected macroscopic cross sections.\n"
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Arguments:\n"
           "    Et    transport corrected total cross section\n"
           "    Ea    absorption cross section\n"
           "    Es    transport corrected scattering cross section matrix",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"))

      .def("ngroups", &TransportXS::ngroups, "Number of energy groups")

      .def("fissile", &TransportXS::fissile, "True if material is fissile")

      .def("Et", py::overload_cast<>(&TransportXS::Et, py::const_),
           "Transport corrected total cross section array")

      .def("Et", py::overload_cast<std::size_t>(&TransportXS::Et, py::const_),
           "Transport corrected total cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Ea", &TransportXS::Ea,
           "Absorption cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("vEf", &TransportXS::vEf,
           "Fission yield * cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Er", &TransportXS::Er,
           "Removal cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("chi", &TransportXS::chi,
           "Fission spectrum in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Es", py::overload_cast<std::size_t>(&TransportXS::Es, py::const_),
           "Transport corrected scattering cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Es",
           py::overload_cast<std::size_t, std::size_t>(&TransportXS::Es,
                                                       py::const_),
           "Transport corrected scattering cross section from group gin to "
           "gout\n\n"
           "Arguments:\n"
           "    gin   incoming energy group\n"
           "    gout  outgoing energy group",
           py::arg("gin"), py::arg("gout"))

      .def("__mul__", &TransportXS::operator*)
      .def("__add__", &TransportXS::operator+)
      .def("__imul__", &TransportXS::operator*=)
      .def("__iadd__", &TransportXS::operator+=);
}
