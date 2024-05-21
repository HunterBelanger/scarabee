#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/material.hpp>
#include <data/nd_library.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_Nuclide(py::module& m) {
  py::class_<Nuclide>(m, "Nuclide")
      .def(py::init<>(),
           "A name-fraction pair, representing a nuclide in a material "
           "composition")
      .def_readwrite("name", &Nuclide::name,
                     "Name of the nuclide (i.e. U235, H1_H2O, etc.)")

      .def_readwrite("fraction", &Nuclide::fraction,
                     "Fraction of the material (by atoms or weight) that is "
                     "occupied by this nuclide");
}

void init_MaterialComposition(py::module& m) {
  py::enum_<Fraction>(m, "Fraction")
      .value("Atoms", Fraction::Atoms)
      .value("Weight", Fraction::Weight);

  py::class_<MaterialComposition>(m, "MaterialComposition")
      .def(py::init<>(), "Creates an empty material composition")

      .def_readwrite("nuclides", &MaterialComposition::nuclides,
                     "List of nuclides in the material as name-fraction pairs")

      .def_readwrite("fractions", &MaterialComposition::fractions,
                     "Flag indicating if the fractions are in atoms or weight")

      .def("add_nuclide",
           py::overload_cast<const std::string&, double>(
               &MaterialComposition::add_nuclide),
           "Adds a new nuclide to the material.\n\n"
           "Arguments:\n"
           "    name     Name of the nuclide\n"
           "    fraction Fraction that the nuclide occupies in the material",
           py::arg("name"), py::arg("fraction"))

      .def("add_nuclide",
           py::overload_cast<const Nuclide&>(&MaterialComposition::add_nuclide),
           "Adds a new nuclide to the material.\n\n"
           "Arguments:\n"
           "    comp MaterialComponent giving the nuclide name and fraction",
           py::arg("comp"));
}

void init_Material(py::module& m) {
  py::enum_<DensityUnits>(m, "DensityUnits")
      .value("g_cm3", DensityUnits::g_cm3, "grams per cubic-centimeter")
      .value("a_bcm", DensityUnits::a_bcm, "atoms per barn-centimeter")
      .value("sum", DensityUnits::sum,
             "computer density from sum of fractions");

  py::class_<Material>(m, "Material")
      .def(py::init<const MaterialComposition&, double,
                    std::shared_ptr<NDLibrary>>(),
           "Creates a new Material definition.\n\n"
           "Arguments:\n"
           "    comp MaterialComposition defining material components\n"
           "    temp Temperature of the material in kelvin\n"
           "    ndl  NDLibrary instance for nuclear data",
           py::arg("comp"), py::arg("temp"), py::arg("ndl"))

      .def(py::init<const MaterialComposition&, double, double, DensityUnits,
                    std::shared_ptr<NDLibrary>>(),
           "Creates a new Material definition.\n\n"
           "Arguments:\n"
           "    comp MaterialComposition defining material components\n"
           "    temp Temperature of the material in kelvin\n"
           "    density Density of the material in units given by du\n"
           "    du      Units of the provided density (if sum, density is "
           "ignored)\n"
           "    ndl  NDLibrary instance for nuclear data",
           py::arg("comp"), py::arg("temp"), py::arg("density"), py::arg("du"),
           py::arg("ndl"))

      .def("has_component", &Material::has_component,
           "Returns True if the indicated nuclide is present in the "
           "material.\n\n"
           "Arguments:\n"
           "    name  Name of the nuclide",
           py::arg("name"))

      .def("atom_density", &Material::atom_density,
           "Returns the number of atoms per barn-centimeter of the indicated "
           "nuclide.\n\n"
           "Arguments:\n"
           "    name  Name of the nuclide",
           py::arg("name"))

      .def("carlvik_xs", &Material::carlvik_xs,
           "Returns the macroscopic material cross section, self-shielded "
           "according to the Carlvik two-term approximation.\n\n"
           "Arguments:\n"
           "    C    Dancoff correction factor\n"
           "    Ee   Escpae cross section\n"
           "    ndl  Nuclear data library",
           py::arg("C"), py::arg("Ee"), py::arg("ndl"))

      .def("roman_xs", &Material::roman_xs,
           "Returns the macroscopic material cross section, self-shielded "
           "according to the Roman two-term approximation.\n\n"
           "Arguments:\n"
           "    C    Dancoff correction factor\n"
           "    Ee   Escpae cross section\n"
           "    ndl  Nuclear data library",
           py::arg("C"), py::arg("Ee"), py::arg("ndl"))

      .def("dilution_xs", &Material::dilution_xs,
           "Returns the macroscopic material cross section with nuclides "
           "interpolated to the provided dilutions.\n\n"
           "Arguments:\n"
           "    dils  List of dilutions\n"
           "    ndl   Nuclear data library",
           py::arg("dils"), py::arg("ndl"))

      .def_property_readonly(
          "composition", &Material::composition,
          "The MaterialComposition defining the nuclides in the material")

      .def_property_readonly("size", &Material::size,
                             "Number of nuclides in the material")

      .def_property_readonly("temperature", &Material::temperature,
                             "Temperature of the material in kelvin")

      .def_property_readonly("average_molar_mass",
                             &Material::average_molar_mass,
                             "Average molar of an atom in the material, mass "
                             "based on all nuclides in the material")

      .def_property_readonly(
          "atoms_per_bcm", &Material::atoms_per_bcm,
          "Total number of atoms per barn-centimeter in the material")

      .def_property_readonly(
          "potential_xs", &Material::potential_xs,
          "Macroscopic potential scattering cross section in units of 1/cm")

      .def_property_readonly(
          "grams_per_cm3", &Material::grams_per_cm3,
          "Density of the materil in grams per cubic-centimeter")

      .def_property_readonly("fissile", &Material::fissile,
                             "True if the material is fissile, False otherwise")

      .def_property_readonly(
          "resonant", &Material::resonant,
          "True if the material is resonant, False otherwise");
}