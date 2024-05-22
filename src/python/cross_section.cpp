#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include <cross_section.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CrossSection(py::module& m) {
  py::class_<CrossSection, std::shared_ptr<CrossSection>>(m, "CrossSection")
      .def(py::init<const xt::xtensor<double, 1>& /*Etr*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es_tr*/,
                    const xt::xtensor<double, 1>& /*Ef*/,
                    const xt::xtensor<double, 1>& /*vEf*/,
                    const xt::xtensor<double, 1>& /*chi*/,
                    const std::string& /*name*/>(),
           "Transport corrected macroscopic cross sections.\n\n"
           "Arguments:\n"
           "    Etr    transport corrected total cross section\n"
           "    Ea     absorption cross section\n"
           "    Es_tr  transport corrected scattering cross section matrix\n"
           "    Ef     fission cross section\n"
           "    vEf    fission yield time cross section\n"
           "    chi    fission spectrum\n"
           "    name   name of material",
           py::arg("Etr"), py::arg("Ea"), py::arg("Es_tr"), py::arg("Ef"),
           py::arg("vEf"), py::arg("chi"), py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*Et*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const xt::xtensor<double, 2>& /*Es1*/,
                    const xt::xtensor<double, 1>& /*Ef*/,
                    const xt::xtensor<double, 1>& /*vEf*/,
                    const xt::xtensor<double, 1>& /*chi*/,
                    const std::string& /*name*/>(),
           "Transport corrected macroscopic cross sections.\n\n"
           "Arguments:\n"
           "    Et    total cross section\n"
           "    Ea    absorption cross section\n"
           "    Es    scattering cross section matrix\n"
           "    Ef    fission cross section\n"
           "    vEf   fission yield time cross section\n"
           "    chi   fission spectrum\n"
           "    name  name of material",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"), py::arg("Es1"),
           py::arg("Ef"), py::arg("vEf"), py::arg("chi"), py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*Etr*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es_tr*/,
                    const std::string& /*name*/>(),
           "Transport corrected macroscopic cross sections.\n"
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Arguments:\n"
           "    Etr    transport corrected total cross section\n"
           "    Ea     absorption cross section\n"
           "    Es_tr  transport corrected scattering cross section matrix\nn"
           "    name   name of material",
           py::arg("Etr"), py::arg("Ea"), py::arg("Es_tr"),
           py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*Et*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const xt::xtensor<double, 2>& /*Es1*/,
                    const std::string& /*name*/>(),
           "Transport corrected macroscopic cross sections.\n"
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Arguments:\n"
           "    Et    total cross section\n"
           "    Ea    absorption cross section\n"
           "    Es    scattering cross section matrix\nn"
           "    name  name of material",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"), py::arg("Es1"),
           py::arg("name") = "")

      .def("ngroups", &CrossSection::ngroups, "Number of energy groups")

      .def_property("name", &CrossSection::name, &CrossSection::set_name, "Name of material")

      .def("fissile", &CrossSection::fissile, "True if material is fissile")

      .def("Etr", py::overload_cast<>(&CrossSection::Etr, py::const_),
           "Transport corrected total cross section array")

      .def("Etr",
           py::overload_cast<std::size_t>(&CrossSection::Etr, py::const_),
           "Transport corrected total cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Et", &CrossSection::Et,
           "Total cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Ea", &CrossSection::Ea,
           "Absorption cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Ef", &CrossSection::Ef,
           "Fission cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("vEf", &CrossSection::vEf,
           "Fission yield * cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("nu", &CrossSection::vEf,
           "Fission yield in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Er", &CrossSection::nu,
           "Removal cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("chi", &CrossSection::chi,
           "Fission spectrum in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Es_tr",
           py::overload_cast<std::size_t>(&CrossSection::Es_tr, py::const_),
           "Transport corrected scattering cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Es_tr",
           py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es_tr,
                                                       py::const_),
           "Transport corrected scattering cross section from group gin to "
           "gout\n\n"
           "Arguments:\n"
           "    gin   incoming energy group\n"
           "    gout  outgoing energy group",
           py::arg("gin"), py::arg("gout"))

      .def(
          "Es1",
          py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es1,
                                                      py::const_),
          "P1 moment of the scattering cross section from group gin to gout\n\n"
          "Arguments:\n"
          "    gin   incoming energy group\n"
          "    gout  outgoing energy group",
          py::arg("gin"), py::arg("gout"))

      .def("Es", py::overload_cast<std::size_t>(&CrossSection::Es, py::const_),
           "Scattering cross section in group g\n\n"
           "Arguments:\n"
           "    g   energy group",
           py::arg("g"))

      .def("Es",
           py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es,
                                                       py::const_),
           "Scattering cross section from group gin to gout\n\n"
           "Arguments:\n"
           "    gin   incoming energy group\n"
           "    gout  outgoing energy group",
           py::arg("gin"), py::arg("gout"))

      .def("__mul__", &CrossSection::operator*)
      .def("__rmul__", [](const CrossSection& xs, double N) { return xs * N; })
      .def("__add__", &CrossSection::operator+)
      .def("__imul__", &CrossSection::operator*=)
      .def("__iadd__", &CrossSection::operator+=);
}
