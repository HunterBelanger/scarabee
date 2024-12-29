#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <assemblies/pins/fuel_pin.hpp>
#include <assemblies/pins/guide_tube.hpp>
#include <assemblies/pins/burnable_poison_pin.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_FuelPin(py::module& m) {
  py::class_<FuelPin, std::shared_ptr<FuelPin>>(
      m, "FuelPin",
      "Represents a fuel pin, contianing its geometric and material "
      "properties.\n\n"
      "Parameters\n"
      "----------\n"
      "fuel : Material\n"
      "    Material representing the fuel composition, temperature, and "
      "density.\n"
      "fuel_radius : float\n"
      "    Radius of the fuel pellet.\n"
      "gap : Material, optional\n"
      "    Material representing the gap composition, temperature, and "
      "density.\n"
      "gap_radius: float, optional\n"
      "    Outer radius of the gap between the fuel pellet and the cladding\n"
      "    (if modeled).\n"
      "clad: Material\n"
      "    Material representing the cladding composition, temperature, and\n"
      "    density.\n"
      "clad_radius : float\n"
      "    Outer radius of the cladding.\n"
      "fuel_rings : int\n"
      "    The number of rings into which the fuel pellet will be divided.\n"
      "    Default value is 1.\n"
      "needs_buffer : bool\n"
      "    Indicates if the fuel pin needs an outer \"buffer\" region to drive "
      "the flux for the spectrum calculation. Default value is False.")
      .def(py::init<std::shared_ptr<Material> /*fuel*/, double /*fuel_radius*/,
                    std::shared_ptr<Material> /*gap*/,
                    std::optional<double> /*gap_radius*/,
                    std::shared_ptr<Material> /*clad*/, double /*clad_radius*/,
                    std::size_t /*fuel_rings*/, bool /*needs_buffer*/>(),
           py::arg("fuel"), py::arg("fuel_radius"), py::arg("gap"),
           py::arg("gap_radius"), py::arg("clad"), py::arg("clad_radius"),
           py::arg("fuel_rings") = 1, py::arg("needs_buffer") = false);
}

void init_GuideTube(py::module& m) {
  py::class_<GuideTube, std::shared_ptr<GuideTube>>(
      m, "GuideTube",
      "Represents an empty guide tube, contianing its geometric and material\n"
      "properties.\n\n"
      "Parameters\n"
      "----------\n"
      "inner_radius : float\n"
      "    Inside radius of the guide tube.\n"
      "outer_radius : float\n"
      "    Outside radius of the guide tube.\n"
      "clad : Material\n"
      "    Material representing the guide tube composition, temperature, and\n"
      "    density (assumed to be the same material as the fuel pin "
      "cladding).\n\n")
      .def(py::init<std::shared_ptr<Material> /*clad*/, double /*inner_radius*/,
                    double /*outer_radius*/>(),
           py::arg("clad"), py::arg("inner_radius"), py::arg("outer_radius"));
}

void init_BurnablePoisonPin(py::module& m) {
  py::class_<BurnablePoisonPin, std::shared_ptr<BurnablePoisonPin>>(
      m, "BurnablePoisonPin",
      "Represents a burnable poison bin which could be of borosilicate glass\n"
      "(BSG), or a wet annular burnable absorber (WABA).\n\n"
      "Parameters\n"
      "----------\n"
      "center: Material\n"
      "    Material at the center of the burnable poison tube.\n"
      "center_radius: float\n"
      "    Radius of the center material in the burnable poison tube.\n"
      "poison_clad: Material\n"
      "    Cladding material for the burnable poison tube.\n"
      "inner_poison_clad_radius: float\n"
      "    Radius of the inner cladding of the burnable poison tube.\n"
      "gap: Material, optional\n"
      "    Material of the gap between the burnable poison and the cladding\n"
      "    (if present).\n"
      "inner_gap_radius: float, optional\n"
      "    Radius of the inner gap between poison cladding and the poison.\n"
      "poison: Material\n"
      "    Material of the burnable poison.\n"
      "poison_radius: float\n"
      "    Outer radius of the burnable poison.\n"
      "outer_gap_radius: float, optional\n"
      "    Outer radius of the gap between the burnable poison and cladding.\n"
      "outer_poison_clad_radius: float\n"
      "    Outer radius of the cladding of the burnable poison tube.\n"
      "inner_moderator_radius: float\n"
      "    The outer radius of the moderator between the burnable poison tube\n"
      "    and the guide tube.\n"
      "guide_tube_clad: Material\n"
      "    Material for the guide tube.\n"
      "guide_tube_radius: float\n"
      "    Outer radius of the guide tube.\n\n")
      .def(py::init<
               std::shared_ptr<Material> /*center*/, double /*center_radius*/,
               std::shared_ptr<Material> /*poison_clad*/,
               double /*inner_poison_clad_radius*/,
               std::shared_ptr<Material> /*gap*/,
               std::optional<double> /*inner_gap_radius*/,
               std::shared_ptr<Material> /*poison*/, double /*poison_radius*/,
               std::optional<double> /*outer_gap_radius*/,
               double /*outer_poison_clad_radius*/,
               double /*inner_moderator_radius*/,
               std::shared_ptr<Material> /*guide_tube_clad*/,
               double /*guide_tube_radius*/>(),
           py::arg("center"), py::arg("center_radius"), py::arg("poison_clad"),
           py::arg("inner_poison_clad_radius"), py::arg("gap"),
           py::arg("inner_gap_radius"), py::arg("poison"),
           py::arg("poison_radius"), py::arg("outer_gap_radius"),
           py::arg("outer_poison_clad_radius"),
           py::arg("inner_moderator_radius"), py::arg("guide_tube_clad"),
           py::arg("guide_tube_radius"));
}