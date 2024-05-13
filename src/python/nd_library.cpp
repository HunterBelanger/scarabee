#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/nd_library.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_NuclideHandle(py::module& m) {
  py::class_<NuclideHandle>(m, "NuclideHandle")
      .def_readonly("name", &NuclideHandle::name, "Indentifier of the nuclide")

      .def_readonly("label", &NuclideHandle::label,
                    "Optional label provided at library creation")

      .def_readonly(
          "temperatures", &NuclideHandle::temperatures,
          "List of temperatures at which cross sections are tabulated")

      .def_readonly("dilutions", &NuclideHandle::dilutions,
                    "List of dilutions at which cross sections are tabulated")

      .def_readonly("awr", &NuclideHandle::awr,
                    "Atomic weight ratio of the nuclide")

      .def_readonly("potential_xs", &NuclideHandle::potential_xs,
                    "Potential scattering cross section of the nuclide")

      .def_readonly("ZA", &NuclideHandle::ZA,
                    "The ZA number of the nuclide, constructed as Z*1000 + A")

      .def_readonly("fissile", &NuclideHandle::fissile,
                    "True if the nuclide is fissile, False otherwise")

      .def_readonly("resonant", &NuclideHandle::resonant,
                    "True if the nuclide is resonant, False otherwise");
}

void init_NDLibrary(py::module& m) {
  py::class_<NDLibrary, std::shared_ptr<NDLibrary>>(m, "NDLibrary")
  .def(py::init<const std::string&>(),
  "Creates a new Nuclear Data Library object.\n\n"
  "Arguments:\n"
  "    fname  Name of the hdf5 file with the library", py::arg("fname"))

  .def("get_nuclide", py::overload_cast<const std::string&>(&NDLibrary::get_nuclide, py::const_),
  "Returns the NuclideHandle of the the indicated nuclide.\n\n"
  "Arguments:\n"
  "    name  Name of the desired nuclide", py::arg("name"))

  .def("interp_nuclide_xs", &NDLibrary::interp_nuclide_xs,
  "Interpolates the cross section of the prescribed nuclide to the\n"
  "desired temperature and dilution.\n\n"
  "Arguments:\n"
  "    name  Name of the desired nuclide\n"
  "    temp  Desired temperature in kelvin\n"
  "    dil   Desired dilution in barns", py::arg("name"), py::arg("temp"), py::arg("dil"))

  .def("carlvik_two_term", &NDLibrary::carlvik_two_term,
  "Uses Carlvik's two-term rational approximation for self shielding of cross sections.\n\n"
  "Arguments:\n"
  "    name        Name of the nuclide to be treated\n"
  "    mat_pot_xs  Macroscopic potential scattering cross section of the material\n"
  "    temp        Temperature of the material (in kelvin)\n"
  "    N           Atom density of the nuclide\n"
  "    C           Dancoff correction factor\n"
  "    Ee          Escape cross section of the fuel lump",
  py::arg("name"), py::arg("mat_pot_xs"), py::arg("temp"), py::arg("N"), py::arg("C"), py::arg("Ee"))

  .def_property_readonly("library", &NDLibrary::library,
  "Name of the nuclear data library (if provided)")

  .def_property_readonly("ngroups", &NDLibrary::ngroups,
  "Number of energy groups in the library")

  .def_property_readonly("group_bounds", &NDLibrary::group_bounds,
  "The boundaries of the energy groups for the group structure (in decreasing order)")

  .def_property_readonly("group_structure", &NDLibrary::group_structure,
  "The name of the group structure (if provided)");
}