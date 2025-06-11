#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <moc/cmfd.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CMFD(py::module& m) {

  py::enum_<CMFDSurfaceCrossing::Type>(m,"CMFDSurfaceCrossingType")
    .value("XN",CMFDSurfaceCrossing::Type::XN)
    .value("XP",CMFDSurfaceCrossing::Type::XP)
    .value("YN",CMFDSurfaceCrossing::Type::YN)
    .value("YP",CMFDSurfaceCrossing::Type::YP)
    .value("TR",CMFDSurfaceCrossing::Type::TR)
    .value("BR",CMFDSurfaceCrossing::Type::BR)
    .value("BL",CMFDSurfaceCrossing::Type::BL)
    .value("TL",CMFDSurfaceCrossing::Type::TL);
  
  py::class_<CMFDSurfaceCrossing>(m, "CMFDSurfaceCrossing")
  .def_readonly("is_valid", &CMFDSurfaceCrossing::is_valid,
  "Bool that indicates whether or not the segment end is on a valid CMFD surface")
  .def_readonly("cell_index", &CMFDSurfaceCrossing::cell_index,
  "Int representing the linear cell index")
  .def_readonly("crossing", &CMFDSurfaceCrossing::crossing,
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
    
  .def_property(
          "keff_tolerance", &CMFD::keff_tolerance,
          &CMFD::set_keff_tolerance,
          "Maximum relative absolute difference in keff for convergence")

  .def_property(
          "flux_tolerance", &CMFD::flux_tolerance,
          &CMFD::set_flux_tolerance,
          "Maximum relative absolute difference in flux for convergence")

  .def_property(
          "flux_limiting", nullptr,
          &CMFD::set_flux_limiting,
          "Whether or not to use the flux-limiting condition when calculating"
          "the surface and non-linear diffusion coefficients")

  .def_property(
          "larsen_correction", nullptr,
          &CMFD::set_larsen_correction,
          "Whether or not to use Larsen's corrected diffusion coefficient"
          "for optically thick meshes. May reduce performance if not needed")

  .def_property(
          "damping", nullptr,
          &CMFD::set_damping,
          "set CMFD damping factor for under-relaxing the nonlinear diffusion coefficient"
          "between iterations.")

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
      py::arg("point"),py::arg("direction"))
      
  .def("get_tile", &CMFD::get_tile,
      "Finds the tile that a segement end exists in\n\n"
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
      py::arg("point"),py::arg("direction"))
  .def("tally_current",&CMFD::tally_current,
      "Tallies the MOC current onto the appropriate CMFD surface(s)\n\n"
      "Parameters\n"
      "----------\n"
      "aflx: float\n"
      "      value of the angular flux to be tallied"
      "u: Direction\n"
          "Direction of the flux\n"
      "G: unsigned int\n"
      "   Int which gives the # of the flux energy group"
      "surf: CMFDSurfaceCrossing"
      "      CMFDSurfaceCrossing object which gives the information about the surface to tally on\n"
      "Returns\n"
      "-------\n"
      "None",
      py::arg("aflx"),py::arg("u"),py::arg("G"),py::arg("surf"))
  .def("current",
      py::overload_cast<const std::size_t, const std::size_t>(
      &CMFD::current, py::const_),
      "TODO: Fill in docstring",
      py::arg("G"),py::arg("surface"))
  .def("current",
      py::overload_cast<const std::size_t, const std::size_t>(
        &CMFD::current),
      "TODO: Fill in docstring",
      py::arg("G"),py::arg("surface"))
   .def("flux",&CMFD::flux,
      "Gets the CMFD flux at cell i,j in group g \n\n"
      "Parameters\n"
      "----------\n"
      "i: int\n"
      "      cell x index\n"
      "j: int\n"
      "      cell y index\n"
      "g: int\n"
      "      cmfd group\n"
      "Returns\n"
      "-------\n"
      "flux: double\n"
      "      CMFD scalar flux at cell i,j in group g\n",
      py::arg("i"),py::arg("j"),py::arg("g"));
}