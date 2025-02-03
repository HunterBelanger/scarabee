#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/cmfd.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CMFD(py::module& m) {

  py::enum_<CMFDSurfaceCrossing::Type>(m,"CMFDSurfaceCrossingType")
    .value("XN",CMFDSurfaceCrossing::Type::XN)
    .value("XN",CMFDSurfaceCrossing::Type::XP)
    .value("XN",CMFDSurfaceCrossing::Type::YN)
    .value("XN",CMFDSurfaceCrossing::Type::YP)
    .value("XN",CMFDSurfaceCrossing::Type::TR)
    .value("XN",CMFDSurfaceCrossing::Type::BR)
    .value("XN",CMFDSurfaceCrossing::Type::BL)
    .value("XN",CMFDSurfaceCrossing::Type::TL);
  
  py::class_<CMFDSurfaceCrossing>(m, "CMFDSurfaceCrossing")
  .def_property_readonly("is_valid", &CMFDSurfaceCrossing::is_valid,
  "Bool that indicates whether or not the segment end is on a valid CMFD surface")
  .def_property_readonly("cell_index", &CMFDSurfaceCrossing::cell_index,
  "Int with the CMFD cell index")
  .def_property_readonly("crossing", &CMFDSurfaceCrossing::crossing,
  "Defines which corner / side of the cell the segment end passes through");

  py::class_<CMFD, std::shared_ptr<CMFD>>(m, "CMFD")
  .def(py::init<
          const std::vector<double>& /*dx*/, const std::vector<double>& /*dy*/,
          const std::vector<std::pair<std::size_t, std::size_t>>& /*groups*/>(),
      "A Cartesian mesh for accelerating MOC convergence.\n\n"
      "Parameters\n"
      "----------\n"
      "dx : list of float\n"
      "     List of all x widths.\n"
      "dy : list of float\n"
      "     List of all y heights.\n"
      "groups : list of 2D tuples of ints.\n"
      "         The scheme for condensing energy groups.\n",
      py::arg("dx"), py::arg("dy"), py::arg("groups"))
  .def("get_surface", &CMFD::get_surface,
      "Assigns the CMFD surface info to a segement end if it exists\n\n"
      "Parameters\n"
      "----------\n"
      "point: list of float\n"
      "       list of x and y position of segment end\n"
      "direction: list of float\n"
      "           list of x and y components of segment direction unit vector\n"
      "Returns\n"
      "-------\n"
      "surface: CMFDSurfaceCrossing object\n"
      "         Contains the surface crossing information for the segement end\n",
      py::arg("point"),py::arg("direction"));
}