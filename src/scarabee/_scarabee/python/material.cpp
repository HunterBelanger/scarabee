#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/material.hpp>
#include <data/nd_library.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_Nuclide(py::module& m) {
  py::class_<Nuclide>(
      m, "Nuclide",
      "A Nuclide represents a nuclide name - fraction pair, which acts as a "
      "component of a :py:class:`MaterialComposition`.")

      .def(py::init<>(),
           "Creates a Nuclide with no name, and a fraction of zero.")

      .def_readwrite("name", &Nuclide::name,
                     "Name of the nuclide (i.e. \"U235\", \"H1_H2O\", etc.).")

      .def_readwrite("fraction", &Nuclide::fraction,
                     "Fraction of the material (by atoms or weight) that is "
                     "occupied by this nuclide.");
}

void init_MaterialComposition(py::module& m) {
  py::enum_<Fraction>(m, "Fraction")
      .value("Atoms", Fraction::Atoms, "Indicates fractions are in Atoms.")
      .value("Weight", Fraction::Weight, "Indicates fractions are in Weight.");

  py::class_<MaterialComposition, std::shared_ptr<MaterialComposition>>(
      m, "MaterialComposition",
      "A MaterialComposition object represents all nuclides which make up a "
      "material, and portion of the material taken up b y each nuclide.")

      .def(py::init<Fraction, const std::string&>(),
           "Creates an empty composition, with no nuclides.\n\n"
           "Parameters\n"
           "----------\n"
           "fractions : Fraction\n"
           "  Indicates if material is specified using atom or weight "
           "fractions. Default value is Fraction.Atom.\n"
           "name : str\n"
           "  Name for the material (default value is \"\").\n\n",
           py::arg("fractions") = Fraction::Atoms, py::arg("name") = "")

      .def_readonly(
          "nuclides", &MaterialComposition::nuclides,
          "List of :py:class:`Nuclide` objects, defining the compositon.")

      .def_readonly(
          "fractions", &MaterialComposition::fractions,
          "Flag indicating if the nuclide fractions are in Atoms or Weight.")

      .def_readwrite("name", &MaterialComposition::name,
                     "String with the name of the material.")

      .def("add_element", &MaterialComposition::add_element,
           "Adds all naturally occurring isotopes of an element to the "
           "material.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the element.\n"
           "fraction : float\n"
           "       Fraction that the element occupies in the material.\n\n",
           py::arg("name"), py::arg("fraction"))

      .def("add_leu", &MaterialComposition::add_leu,
           "Adds all naturally occurring isotopes of Uranium to the material "
           "for a given weight percent enrichment. Method is only valid for "
           "enrichments <= 5 w/o. More information can be found in "
           "`ORNL/CSD/TM-244 <https://doi.org/10.2172/5561567>`_.\n\n"
           "Parameters\n"
           "----------\n"
           "enrichment : float\n"
           "       Weight enrichment of U235.\n"
           "fraction : float\n"
           "       Fraction that the Uranium occupies in the material.\n\n",
           py::arg("enrichment"), py::arg("fraction"))

      .def("add_nuclide",
           py::overload_cast<const std::string&, double>(
               &MaterialComposition::add_nuclide),
           "Adds a new nuclide to the material.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the nuclide.\n"
           "fraction : float\n"
           "       Fraction that the nuclide occupies in the material.\n\n",
           py::arg("name"), py::arg("fraction"))

      .def(
          "add_nuclide",
          py::overload_cast<const Nuclide&>(&MaterialComposition::add_nuclide),
          "Adds a new nuclide to the material.\n\n"
          "Parameters\n"
          "----------\n"
          "nuc : Nuclide\n"
          "      :py:class:`Nuclide` giving the nuclide name and fraction.\n\n",
          py::arg("nuc"))

      .def("__deepcopy__", [](const MaterialComposition& comp) {
        return MaterialComposition(comp);
      });
}

void init_Material(py::module& m) {
  py::enum_<DensityUnits>(m, "DensityUnits")
      .value("g_cm3", DensityUnits::g_cm3, "Grams per cubic-centimeter.")
      .value("a_bcm", DensityUnits::a_bcm, "Atoms per barn-centimeter.")
      .value("sum", DensityUnits::sum,
             "Compute density from sum of fractions.");

  py::class_<Material, std::shared_ptr<Material>>(
      m, "Material",
      "A Material object represents a material composition in combination with "
      "a density and temperature. From a Material, it is possible to obtain "
      "macroscopic potential cross sections, and construct self-shielded "
      "macroscopic cross sections.")

      .def(py::init<const MaterialComposition&, double,
                    std::shared_ptr<NDLibrary>>(),
           "Creates a new Material definition with specified temperature. The "
           "density is determined by computing the sum of fractions in the "
           ":py:class:`MaterialCompositon`. If the fractions are atoms, then "
           "the density units are interpreted as atoms per barn-centimeter. "
           "If the fractions are weight, then the density units are "
           "interpreted as grams per cubic-centimeter.\n\n"
           "Parameters\n"
           "----------\n"
           "comp : MaterialComposition\n"
           "       Composition defining the material.\n"
           "temp : float\n"
           "       Temperature of the material in kelvin.\n"
           "ndl  : NDLibrary\n"
           "       Nuclear data library.\n\n",
           py::arg("comp"), py::arg("temp"), py::arg("ndl"))

      .def(py::init<const MaterialComposition&, double, double, DensityUnits,
                    std::shared_ptr<NDLibrary>>(),
           "Creates a new Material definition with specified temperature and "
           "density. If the density units are sum, then the provided value of "
           "the density is ignored, and instead computed from the composition "
           "fractions.\n\n"
           "Parameters\n"
           "----------\n"
           "comp : MaterialComposition\n"
           "       Composition defining the material.\n"
           "temp : float\n"
           "       Temperature of the material in kelvin.\n"
           "density : float\n"
           "          Density of the material in units given by du.\n"
           "du : DensityUnits\n"
           "     Units of the provided density (if sum, density is ignored).\n"
           "ndl : NDLibrary\n"
           "      Nuclear data library.\n\n",
           py::arg("comp"), py::arg("temp"), py::arg("density"), py::arg("du"),
           py::arg("ndl"))

      .def("has_component", &Material::has_component,
           "Indicates if a nuclide is present in the material.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the desired nuclide.\n\n"
           "Returns\n"
           "-------\n"
           "bool\n"
           "     True if nuclide is present, False if not.",
           py::arg("name"))

      .def("atom_density", &Material::atom_density,
           "Number of atoms per barn-centimeter of the indicated nuclide.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the desired nuclide.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "      Number of atoms per barn-centimeter of desired nuclide.",
           py::arg("name"))

      .def("carlvik_xs", &Material::carlvik_xs,
           "Computes the macroscopic material cross section, self-shielded "
           "according to the Carlvik two-term approximation.\n\n"
           "Parameters\n"
           "----------\n"
           "C : float\n"
           "    Dancoff correction factor.\n"
           "Ee : float\n"
           "     Escpae cross section.\n"
           "ndl : NDLibrary\n"
           "      Nuclear data library for cross section interpolation.\n"
           "max_l : int\n"
           "        Maximum legendre moment (default is 1).\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "             The macroscopic self-shielded cross section.",
           py::arg("C"), py::arg("Ee"), py::arg("ndl"), py::arg("max_l") = 1)

      .def("roman_xs", &Material::roman_xs,
           "Computes the macroscopic material cross section, self-shielded "
           "according to the Roman two-term approximation.\n\n"
           "Parameters\n"
           "----------\n"
           "C : float\n"
           "    Dancoff correction factor.\n"
           "Ee : float\n"
           "     Escpae cross section.\n"
           "ndl : NDLibrary\n"
           "      Nuclear data library for cross section interpolation.\n"
           "max_l : int\n"
           "        Maximum legendre moment (default is 1).\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "             The macroscopic self-shielded cross section.",
           py::arg("C"), py::arg("Ee"), py::arg("ndl"), py::arg("max_l") = 1)

      .def("dilution_xs", &Material::dilution_xs,
           "Computes the macroscopic material cross section with nuclides "
           "interpolated to the provided dilutions.\n\n"
           "Parameters\n"
           "----------\n"
           "dils : list of float\n"
           "       Desired dilution for each nuclide.\n"
           "ndl : NDLibrary\n"
           "      Nuclear data library for cross section interpolation.\n"
           "max_l : int\n"
           "        Maximum legendre moment (default is 1).\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "             The macroscopic cross section.",
           py::arg("dils"), py::arg("ndl"), py::arg("max_l") = 1)

      .def("ring_carlvik_xs", &Material::ring_carlvik_xs,
           "Computes the macroscopic material cross section, self-shielded "
           "according to the Carlvik two-term approximation for a single ring "
           "of fuel using the Stoker-Weiss method.\n\n"
           "Parameters\n"
           "----------\n"
           "C : float\n"
           "    Dancoff correction factor.\n"
           "Rfuel : float\n"
           "     Radius of the fuel pellet.\n"
           "Rin : float\n"
           "     Inner radius of the fuel ring.\n"
           "Rout : float\n"
           "     Outer radius of the fuel ring.\n"
           "ndl : NDLibrary\n"
           "      Nuclear data library for cross section interpolation.\n"
           "max_l : int\n"
           "        Maximum legendre moment (default is 1).\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "             The macroscopic self-shielded cross section.",
           py::arg("C"), py::arg("Rfuel"), py::arg("Rin"), py::arg("Rout"),
           py::arg("ndl"), py::arg("max_l") = 1)

      .def("clear_micro_xs_data", &Material::clear_micro_xs_data,
           "Clears all of the previously computed microscopic cross section "
           "data.")

      .def_property("max_legendre_order", &Material::max_legendre_order,
                    &Material::set_max_legendre_order,
                    "The maximum legendre order for loading and interpolating "
                    "scattering matrices. The default value is 1.")

      .def_property_readonly("composition", &Material::composition,
                             "The :py:class:`MaterialComposition` defining the "
                             "nuclides in the material.")

      .def_property_readonly("size", &Material::size,
                             "Number of nuclides in the material.")

      .def_property("temperature", &Material::temperature,
                    &Material::set_temperature,
                    "Temperature of the material in kelvin.")

      .def_property_readonly("average_molar_mass",
                             &Material::average_molar_mass,
                             "Average molar of an atom in the material, "
                             "based on all nuclides in the material.")

      .def_property_readonly(
          "atoms_per_bcm", &Material::atoms_per_bcm,
          "Total number of atoms per barn-centimeter in the material.")

      .def_property_readonly(
          "potential_xs", &Material::potential_xs,
          "Macroscopic potential scattering cross section in units of 1/cm.")

      .def_property_readonly(
          "grams_per_cm3", &Material::grams_per_cm3,
          "Density of the materil in grams per cubic-centimeter.")

      .def_property_readonly(
          "fissile", &Material::fissile,
          "True if the material is fissile, False otherwise.")

      .def_property_readonly(
          "resonant", &Material::resonant,
          "True if the material is resonant, False otherwise.")

      .def_property("name", &Material::name, &Material::set_name,
                    "String with the name of the Material.")

      .def("__deepcopy__", [](const Material& mat) { return Material(mat); });

  py::enum_<MixingFraction>(m, "MixingFraction")
      .value("Atoms", MixingFraction::Atoms,
             "Indicates fractions are in Atoms.")
      .value("Weight", MixingFraction::Weight,
             "Indicates fractions are in Weight.")
      .value("Volume", MixingFraction::Weight,
             "Indicates fractions are in Volume.");

  m.def("mix_materials", &mix_materials,
        "Creates a new material defined as a mixture of materials.\n\n"
        "Parameters\n"
        "----------\n"
        "mats : list of Material\n"
        "       Materials in the mixture.\n"
        "fracs : list of float\n"
        "       Fraction for each material.\n"
        "f : MixingFraction\n"
        "    Indicates if provided fractions are in Atoms, Weight, or Volume.\n"
        "ndl : NDLibrary\n"
        "      Nuclear data library.\n\n"
        "Returns\n"
        "-------\n"
        "Material\n"
        "  Mixture material with averaged temperature.\n",
        py::arg("mats"), py::arg("fracs"), py::arg("f"), py::arg("ndl"));
}
