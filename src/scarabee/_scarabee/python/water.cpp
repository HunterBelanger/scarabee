#include <pybind11/pybind11.h>

#include <data/water.hpp>
#include <data/nd_library.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_WaterFuncs(py::module& m) {
  m.def(
      "water_density", &water_density,
      "Computes the density of water at a desired temperature and pressure.\n\n"
      "Parameters\n"
      "----------\n"
      "temperature : float\n"
      "  Temperature of the water in Kelvin.\n"
      "pressure : float\n"
      "  Pressure of the water in MPa.\n\n"
      "Returns\n"
      "-------\n"
      "float\n"
      "  Density of the water in g/cm^3.\n",
      py::arg("temperature"), py::arg("pressure"));

  m.def("borated_water", &borated_water,
        "Makes a :py:class:`Material` for water at a desired temperature, "
        "pressure, and boron concentration.\n\n"
        "Parameters\n"
        "----------\n"
        "boron_ppm : float\n"
        "  Concentration of boron in parts per million.\n"
        "temperature : float\n"
        "  Temperature of the water in Kelvin.\n"
        "pressure : float\n"
        "  Pressure of the water in MPa.\n"
        "ndl : NDLibrary\n"
        "  Nuclear data library.\n\n"
        "Returns\n"
        "-------\n"
        "Material\n"
        "  Material which contains borated water at the desired temperature, "
        "pressure, and concentration.\n",
        py::arg("boron_ppm"), py::arg("temperature"), py::arg("pressure"),
        py::arg("ndl"));
}