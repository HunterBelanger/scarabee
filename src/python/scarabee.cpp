#include <pybind11/pybind11.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>

#include <utils/version.hpp>

namespace py = pybind11;

extern void init_CrossSection(py::module&);
extern void init_DiffusionCrossSection(py::module&);
extern void init_NuclideHandle(py::module&);
extern void init_NDLibrary(py::module&);
extern void init_Nuclide(py::module&);
extern void init_MaterialComposition(py::module&);
extern void init_Material(py::module&);
extern void init_CylindricalCell(py::module&);
extern void init_CylindricalFluxSolver(py::module&);
extern void init_Logging(py::module&);
extern void init_Vector(py::module&);
extern void init_Direction(py::module&);
extern void init_PolarQuadrature(py::module&);
extern void init_BoundaryCondition(py::module&);
extern void init_SimulationMode(py::module&);
extern void init_Track(py::module&);
extern void init_Cell(py::module&);
extern void init_EmptyCell(py::module&);
extern void init_SimplePinCell(py::module&);
extern void init_PinCell(py::module&);
extern void init_Cartesian2D(py::module&);
extern void init_MOCDriver(py::module&);
extern void init_CriticalitySpectrum(py::module&);
extern void init_DiffusionGeometry(py::module& m);

PYBIND11_MODULE(scarabee, m) {
  xt::import_numpy();

  init_Logging(m);
  init_Vector(m);
  init_Direction(m);
  init_CrossSection(m);
  init_DiffusionCrossSection(m);
  init_NuclideHandle(m);
  init_NDLibrary(m);
  init_Nuclide(m);
  init_MaterialComposition(m);
  init_Material(m);
  init_CylindricalCell(m);
  init_CylindricalFluxSolver(m);
  init_PolarQuadrature(m);
  init_BoundaryCondition(m);
  init_SimulationMode(m);
  init_Track(m);
  init_Cell(m);
  init_EmptyCell(m);
  init_SimplePinCell(m);
  init_PinCell(m);
  init_Cartesian2D(m);
  init_MOCDriver(m);
  init_CriticalitySpectrum(m);
  init_DiffusionGeometry(m);

  m.attr("__author__") = "Hunter Belanger";
  m.attr("__copyright__") = "Copyright 2024, Hunter Belanger";
  m.attr("__license__") = "GPL-3.0-or-later";
  m.attr("__maintainer__") = "Hunter Belanger";
  m.attr("__email__") = "hunter.belanger@gmail.com";
  m.attr("__version__") = scarabee::VERSION_STRING;
}
