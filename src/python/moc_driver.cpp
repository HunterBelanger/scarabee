#include <pybind11/pybind11.h>

#include <moc/moc_driver.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_MOCDriver(py::module& m) {
  py::class_<MOCDriver>(m, "MOCDriver")
  .def(py::init<std::shared_ptr<Cartesian2D> /*geometry*/,
                BoundaryCondition /*xmin = BoundaryCondition::Reflective*/,
                BoundaryCondition /*xmax = BoundaryCondition::Reflective*/,
                BoundaryCondition /*ymin = BoundaryCondition::Reflective*/,
                BoundaryCondition /*ymax = BoundaryCondition::Reflective*/>(),
  "Initializes a Method of Characteristics problem.\n\n"
  "Arguments:\n"
  "    geometry  a Cartesian2D instance for the problem\n"
  "    xminbc    BoundaryCondition for minimmum x value\n"
  "    xmaxbc    BoundaryCondition for maximum x value\n"
  "    yminbc    BoundaryCondition for minimmum y value\n"
  "    ymaxbc    BoundaryCondition for maximum y value\n",
  py::arg("geometry"), py::arg("xminbc") = BoundaryCondition::Reflective,
  py::arg("xmaxbc") = BoundaryCondition::Reflective, py::arg("yminbc") = BoundaryCondition::Reflective, py::arg("ymaxbc") = BoundaryCondition::Reflective)

  .def("generate_tracks", py::overload_cast<std::uint32_t, double, PolarQuadrature, bool>(&MOCDriver::generate_tracks),
  "Traces tracks for the calculation across the geometry.\n\n"
  "Arguments:\n"
  "    nangles       number of azimuthal angles (even)\n"
  "    d             max spacing between tracks of a given angle\n"
  "    polar_quad    polar quadrature for generating segment lengths\n"
  "    precalc_exps  precalculate exponentials (default = True)",
  py::arg("nangles"), py::arg("d"), py::arg("polar_quad"), py::arg("precalc_exps")=true)
  
  .def("drawn", &MOCDriver::drawn, "return True if geometry has been traced")
  
  .def_property("keff_tolerance", &MOCDriver::keff_tolerance, &MOCDriver::set_keff_tolerance)
  
  .def_property("flux_tolerance", &MOCDriver::flux_tolerance, &MOCDriver::set_flux_tolerance)
  
  .def("keff", &MOCDriver::keff)
  
  .def("ngroups", &MOCDriver::ngroups)

  .def("polar_quadrature", &MOCDriver::polar_quadrature)
  
  .def("solve_keff", &MOCDriver::solve_keff)
  
  .def_property("x_min_bc", [](const MOCDriver& md) -> BoundaryCondition {return md.x_min_bc();}, [](MOCDriver& md) -> BoundaryCondition& {return md.x_min_bc();})

  .def_property("x_max_bc", [](const MOCDriver& md) -> BoundaryCondition {return md.x_max_bc();}, [](MOCDriver& md) -> BoundaryCondition& {return md.x_max_bc();})

  .def_property("y_min_bc", [](const MOCDriver& md) -> BoundaryCondition {return md.y_min_bc();}, [](MOCDriver& md) -> BoundaryCondition& {return md.y_min_bc();})

  .def_property("y_max_bc", [](const MOCDriver& md) -> BoundaryCondition {return md.y_max_bc();}, [](MOCDriver& md) -> BoundaryCondition& {return md.y_max_bc();});


}
