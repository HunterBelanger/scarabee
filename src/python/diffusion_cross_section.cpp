#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include <diffusion_cross_section.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_DiffusionCrossSection(py::module& m) {
  py::class_<DiffusionCrossSection, std::shared_ptr<DiffusionCrossSection>>(
      m, "DiffusionCrossSection",
      "A DiffusionCrossSection object contains all necessary cross section "
      "information to perform a diffusion calculation.")

      .def(py::init<const xt::xtensor<double, 1>& /*D*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const xt::xtensor<double, 1>& /*Ef*/,
                    const xt::xtensor<double, 1>& /*vEf*/,
                    const xt::xtensor<double, 1>& /*chi*/,
                    const std::string& /*name*/>(),
           "Creates a DiffusionCrossSection from transport corrected data.\n\n"
           "Parameters\n"
           "----------\n"
           "D : ndarray\n"
           "      Diffusion coefficients.\n"
           "Ea : ndarray\n"
           "     Absorption cross section.\n"
           "Es : ndarray\n"
           "     Transport corrected scattering cross section matrix.\n"
           "Ef : ndarray\n"
           "     Fission cross section.\n"
           "vEf : ndarray\n"
           "      Fission yield time cross section.\n"
           "chi : ndarray\n"
           "      Fission spectrum.\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
           py::arg("D"), py::arg("Ea"), py::arg("Es"), py::arg("Ef"),
           py::arg("vEf"), py::arg("chi"), py::arg("name") = "")

      .def(py::init<const xt::xtensor<double, 1>& /*D*/,
                    const xt::xtensor<double, 1>& /*Ea*/,
                    const xt::xtensor<double, 2>& /*Es*/,
                    const std::string& /*name*/>(),
           "Creates a DiffusionCrossSection from transport corrected data. "
           "Fission cross sections and spectrum initialized to zero.\n\n"
           "Parameters\n"
           "----------\n"
           "D : ndarray\n"
           "    Diffusion coefficients.\n"
           "Ea : ndarray\n"
           "     Absorption cross section.\n"
           "Es : ndarray\n"
           "     Transport corrected scattering cross section matrix.\n"
           "name : str (optional)\n"
           "       Name of material.\n\n",
           py::arg("D"), py::arg("Ea"), py::arg("Es"), py::arg("name") = "")

      .def_property_readonly("ngroups", &DiffusionCrossSection::ngroups, "Number of energy groups.")

      .def_property("name", &DiffusionCrossSection::name, &DiffusionCrossSection::set_name,
                    "Name of material.")

      .def_property_readonly("fissile", &DiffusionCrossSection::fissile, "True if material is fissile.")

      .def("D", &DiffusionCrossSection::D,
           "Diffusion coefficient in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Ea", &DiffusionCrossSection::Ea,
           "Absorption cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Ef", &DiffusionCrossSection::Ef,
           "Fission cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("vEf", &DiffusionCrossSection::vEf,
           "Fission yield * cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("nu", &DiffusionCrossSection::nu,
           "Fission yield in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Er", &DiffusionCrossSection::Er,
           "Removal cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("chi", &DiffusionCrossSection::chi,
           "Fission spectrum in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Es", py::overload_cast<std::size_t>(&DiffusionCrossSection::Es, py::const_),
           "Transport corrected scattering cross section in group g.\n\n"
           "Parameters\n"
           "----------\n"
           "g : int\n"
           "    Energy group.",
           py::arg("g"))

      .def("Es", py::overload_cast<std::size_t, std::size_t>(&DiffusionCrossSection::Es, py::const_),
           "Transport corrected scattering cross section from group gin to "
           "gout\n\n"
           "Parameters\n"
           "----------\n"
           "gin : int\n"
           "      Incoming energy group.\n"
           "gout : int\n"
           "       Outgoing energy group.",
           py::arg("gin"), py::arg("gout"));
}
