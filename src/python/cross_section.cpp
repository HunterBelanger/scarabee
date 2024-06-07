#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <cross_section.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CrossSection(py::module& m) {
  py::class_<CrossSection, std::shared_ptr<CrossSection>>(
      m, "CrossSection",
      "A CrossSection object contains all necessary cross section information "
      "to perform a transport calculation. A reduced set of arithmetic "
      "operations can be performed on CrossSections to facilitate the "
      "construction of macroscopic material cross sections. This includes "
      "scalar multiplication, and the addition of CrossSection objects.")

      .def(py::init<const xt::xtensor<double, 1>& /*Etr*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es_tr*/,
                    const xt::xtensor<double, 1>& /*Ef*/,
                    const xt::xtensor<double, 1>& /*vEf*/,
                    const xt::xtensor<double, 1>& /*chi*/,
                    const std::string& /*name*/>(),
           "Creates a CrossSection from transport corrected data.\n\n"
           "Parameters\n"
           "----------\n"
           "Etr : ndarray\n"
           "      Transport corrected total cross section.\n"
           "Ea : ndarray\n"
           "     Absorption cross section.\n"
           "Es_tr : ndarray\n"
           "        Transport corrected scattering cross section matrix.\n"
           "Ef : ndarray\n"
           "     Fission cross section.\n"
           "vEf : ndarray\n"
           "      Fission yield time cross section.\n"
           "chi : ndarray\n"
           "      Fission spectrum.\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
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
           "Creates a CrossSection object from uncorrected cross section "
           "data.\n\n"
           "Parameters\n"
           "----------\n"
           "Et : ndarray\n"
           "     Total cross section\n"
           "Ea : ndarray\n"
           "     Absorption cross section\n"
           "Es : ndarray\n"
           "     Scattering cross section matrix\n"
           "Es1 : ndarray\n"
           "     P1 scattering cross section matrix\n"
           "Ef : ndarray\n"
           "     Fission cross section\n"
           "vEf : ndarray\n"
           "      Fission yield time cross section\n"
           "chi : ndarray\n"
           "      Fission spectrum\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"), py::arg("Es1"),
           py::arg("Ef"), py::arg("vEf"), py::arg("chi"), py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*Etr*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es_tr*/,
                    const std::string& /*name*/>(),
           "Creates a CrossSection from transport corrected data. "
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Parameters\n"
           "----------\n"
           "Etr : ndarray\n"
           "      Transport corrected total cross section.\n"
           "Ea : ndarray\n"
           "     Absorption cross section.\n"
           "Es_tr : ndarray\n"
           "        Transport corrected scattering cross section matrix.\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
           py::arg("Etr"), py::arg("Ea"), py::arg("Es_tr"),
           py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*Et*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const xt::xtensor<double, 2>& /*Es1*/,
                    const std::string& /*name*/>(),
           "Creates a CrossSection object from uncorrected cross section data. "
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Parameters\n"
           "----------\n"
           "Et : ndarray\n"
           "     Total cross section.\n"
           "Ea : ndarray\n"
           "     Absorption cross section.\n"
           "Es : ndarray\n"
           "     Scattering cross section matrix\n"
           "Es1 : ndarray\n"
           "      P1 scattering cross section matrix\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
           py::arg("Et"), py::arg("Ea"), py::arg("Es"), py::arg("Es1"),
           py::arg("name") = "")

      .def_property_readonly("ngroups", &CrossSection::ngroups,
                             "Number of energy groups.")

      .def_property("name", &CrossSection::name, &CrossSection::set_name,
                    "Name of material.")

      .def_property_readonly("fissile", &CrossSection::fissile,
                             "True if material is fissile.")

      .def_property_readonly("anisotropic", &CrossSection::anisotropic,
                             "True if the material has a P1 scattering matrix.")

      .def("Etr", py::overload_cast<>(&CrossSection::Etr, py::const_),
           "Transport corrected total cross section array.")

      .def("Etr",
           py::overload_cast<std::size_t>(&CrossSection::Etr, py::const_),
           "Transport corrected total cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Et", &CrossSection::Et,
           "Total cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Ea", &CrossSection::Ea,
           "Absorption cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Ef", &CrossSection::Ef,
           "Fission cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("vEf", &CrossSection::vEf,
           "Fission yield * cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("nu", &CrossSection::vEf,
           "Fission yield in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Er", &CrossSection::nu,
           "Removal cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("chi", &CrossSection::chi,
           "Fission spectrum in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Es_tr",
           py::overload_cast<std::size_t>(&CrossSection::Es_tr, py::const_),
           "Transport corrected scattering cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Es_tr",
           py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es_tr,
                                                       py::const_),
           "Transport corrected scattering cross section from group gin to "
           "gout\n\n"
           "Parameters\n"
           "----------\n"
           "gin : int\n"
           "      Incoming energy group.\n"
           "gout : int\n"
           "       Outgoing energy group.",
           py::arg("gin"), py::arg("gout"))

      .def(
          "Es1",
          py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es1,
                                                      py::const_),
          "P1 moment of the scattering cross section from group gin to gout\n\n"
          "Parameters\n"
          "----------\n"
          "gin : int\n"
          "      Incoming energy group.\n"
          "gout : int\n"
          "       Outgoing energy group.",
          py::arg("gin"), py::arg("gout"))

      .def("Es", py::overload_cast<std::size_t>(&CrossSection::Es, py::const_),
           "Scattering cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Es",
           py::overload_cast<std::size_t, std::size_t>(&CrossSection::Es,
                                                       py::const_),
           "Scattering cross section from group gin to gout.\n\n"
           "Parameters\n"
           "----------\n"
           "gin : int\n"
           "      Incoming energy group.\n"
           "gout : int\n"
           "       Outgoing energy group.",
           py::arg("gin"), py::arg("gout"))

      .def(
          "condense", &CrossSection::condense,
          "Condenses the cross sections to a new energy group structure. The "
          "condensation group structure is provided as a list of pairs "
          "(2D tuples), indicating the lower and upper bounds (inclusive) of "
          "a macro energy group. \n\n"
          "Parameters\n"
          "----------\n"
          "groups : list of 2D tuples of ints.\n"
          "         The scheme for condensing energy groups.\n"
          "flux : list of floats.\n"
          "       The weighting flux spectrum of the original group structure.",
          py::arg("groups"), py::arg("flux"))

      .def("__mul__", &CrossSection::operator*)
      .def("__rmul__", [](const CrossSection& xs, double N) { return xs * N; })
      .def("__add__", &CrossSection::operator+)
      .def("__imul__", &CrossSection::operator*=)
      .def("__iadd__", &CrossSection::operator+=);
}
