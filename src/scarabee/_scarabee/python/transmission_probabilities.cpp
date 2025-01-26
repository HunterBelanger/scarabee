#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <transmission_probabilities.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_TransmissionProbabilities(py::module& m) {
  py::class_<TransmissionProbabilities,
             std::shared_ptr<TransmissionProbabilities>>(
      m, "TransmissionProbabilities")
      .def(
          py::init<
              const std::vector<std::shared_ptr<CrossSection>>& /*xs*/,
              const std::vector<double>& /*dx*/,
              const std::vector<double>& /*dy*/, BoundaryCondition /*x_min_bc*/,
              BoundaryCondition /*x_max_bc*/, BoundaryCondition /*y_min_bc*/,
              BoundaryCondition /*y_max_bc*/>(),
          "Initializesa transmission probabilities problem.\n\n"
          "Parameters\n"
          "----------\n"
          "xs : list of CrossSection\n"
          "     Cross sections for all regions of the problem.\n"
          "dx : list of float\n"
          "     Widths of all regions along x axis.\n"
          "dy : list of float\n"
          "     Height of all regions along y axis.\n"
          "x_min_bc : BoundaryCondition\n"
          "           Boundary condition at the lower x boundary.\n"
          "x_max_bc : BoundaryCondition\n"
          "           Boundary condition at the upper x boundary.\n"
          "y_min_bc : BoundaryCondition\n"
          "           Boundary condition at the lower y boundary.\n"
          "y_max_bc : BoundaryCondition\n"
          "           Boundary condition at the upper y boundary.\n")

      .def_property_readonly("ngroups", &TransmissionProbabilities::ngroups,
                             "Number of energy groups.")

      .def_property_readonly("nregions", &TransmissionProbabilities::nregions,
                             "Number of regions.")

      .def_property_readonly("nx", &TransmissionProbabilities::nx,
                             "Number of regions along the x axis.")

      .def_property_readonly("ny", &TransmissionProbabilities::ny,
                             "Number of regions along the y axis.")

      .def_property_readonly("solved", &TransmissionProbabilities::solved,
                             "True if solve has been run sucessfully (reset to "
                             "False on generate_tracks).")

      .def_property_readonly("keff", &TransmissionProbabilities::keff,
                             "Value of keff estimated by solver (1 by default "
                             "if no solution has been obtained).")

      .def_property(
          "keff_tolerance", &TransmissionProbabilities::keff_tolerance,
          &TransmissionProbabilities::set_keff_tolerance,
          "Maximum relative absolute difference in keff for convergence.")

      .def_property(
          "flux_tolerance", &TransmissionProbabilities::flux_tolerance,
          &TransmissionProbabilities::set_flux_tolerance,
          "Maximum relative absolute difference in flux for convergence.")

      .def("generate_tracks", &TransmissionProbabilities::generate_tracks,
           "Traces tracks across all tiles to compute the transmission "
           "probabilities.\n\n"
           "Parameters\n"
           "----------\n"
           "nangles : int\n"
           "          Number of azimuthal angles (must be even).\n"
           "d : float\n"
           "    Max spacing between tracks of a given angle (in cm).\n",
           py::arg("nangles"), py::arg("d"))

      .def("solve", &TransmissionProbabilities::solve,
           "Begins iterations to solve problem.")

      .def("flux", &TransmissionProbabilities::flux,
           "Returns the scalar flux in group g at tile index (i,j).\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "     Tile index along x.\n"
           "j : int\n"
           "     Tile index along y.\n"
           "g : int\n"
           "     Energy group index.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "     Flux at (i,j) in group g.\n",
           py::arg("i"), py::arg("j"), py::arg("g"))

      .def("flux_spectrum", &TransmissionProbabilities::flux_spectrum,
           "Returns the flux energy spectrum at tile index (i,j).\n\n"
           "Parameters\n"
           "----------\n"
           "i : int\n"
           "     Tile index along x.\n"
           "j : int\n"
           "     Tile index along y.\n\n"
           "Returns\n"
           "-------\n"
           "ndarray of floats\n"
           "     Flux spectrum at (i,j).\n",
           py::arg("i"), py::arg("j"));
}