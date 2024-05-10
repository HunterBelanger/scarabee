#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/material.hpp>
#include <data/nd_library.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_MaterialComponent(py::module& m) {
  py::class_<MaterialComponent>(m, "MaterialComponent")
    .def(py::init<>())
    .def_readwrite("name", &MaterialComponent::name)
    .def_readwrite("fraction", &MaterialComponent::fraction);
}

void init_MaterialComposition(py::module& m) {
  py::enum_<Fraction>(m, "Fraction")
  .value("Atoms", Fraction::Atoms)
  .value("Weight", Fraction::Weight);

  py::class_<MaterialComposition>(m, "MaterialComposition")
  .def(py::init<>())
  .def_readwrite("components", &MaterialComposition::components)
  .def_readwrite("fractions", &MaterialComposition::fractions)
  .def("add_nuclide", py::overload_cast<const std::string&, double>(&MaterialComposition::add_nuclide))
  .def("add_nuclide", py::overload_cast<const MaterialComponent&>(&MaterialComposition::add_nuclide));
}

void init_Material(py::module& m) {
  py::enum_<DensityUnits>(m, "DensityUnits")
  .value("g_cm3", DensityUnits::g_cm3)
  .value("a_bcm", DensityUnits::a_bcm)
  .value("sum", DensityUnits::sum);

  py::class_<Material>(m, "Material")
  .def(py::init<const MaterialComposition&, double, double, DensityUnits, std::shared_ptr<NDLibrary>>())
  .def_property_readonly("composition", &Material::composition)
  .def_property_readonly("temperature", &Material::temperature)
  .def_property_readonly("average_molar_mass", &Material::average_molar_mass)
  .def_property_readonly("atoms_per_bcm", &Material::atoms_per_bcm)
  .def_property_readonly("potential_xs", &Material::potential_xs)
  .def_property_readonly("grams_per_cm3", &Material::grams_per_cm3)
  .def_property_readonly("fissile", &Material::fissile)
  .def_property_readonly("resonant", &Material::resonant)
  .def("has_component", &Material::has_component)
  .def("atom_density", &Material::atom_density);
}