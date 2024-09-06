#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <assemblies/pwr_reflector.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_PWRReflector(py::module& m) {
  py::class_<PWRReflector, std::shared_ptr<PWRReflector>>(
      m, "PWRReflector",
      "A PWRReflector instance is responsible for performing all the lattice\n"
      "calculations necessary to produce few-group cross sections for a \n"
      "PWR reflector, with a gap, baffle, and water reflector.\n\n"
      "Parameters\n"
      "----------\n"
      "pitch : float\n"
      "    Regular spacing between fuel pins.\n"
      "moderator : Material\n"
      "    Material representing the moderator composition, temperature, and\n"
      "    density. This is used for all of the moderator in the assembly.\n"
      "shape : (int, int)\n"
      "    A tuple containing the number of pins across and down for the fuel\n"
      "    assembly.\n"
      "gap_width : float\n"
      "    Width of the water gap between the assembly and the baffle.\n"
      "baffle_width : float\n"
      "    Width of the core baffle. Can be zero if baffle is None.\n"
      "baffle : Material\n"
      "    Material respresenting the core baffle. Can be None is baffle_width "
      "is zero.\n"
      "ndl : NDLibrary\n"
      "    The nuclear data library to be used in the calculations.\n\n"
      "Attributes\n"
      "----------\n"
      "condensation_scheme : list of pairs of ints\n"
      "    Defines how the energy groups will be condensed from the "
      "microgroup\n"
      "    structure of the nuclear data library, to the macrogroup structure\n"
      "    used in the full assembly calculation.\n"
      "few_group_condensation_scheme : list of pairs of ints, optional\n"
      "    Defines how the energy groups will be condensed from the "
      "macrogroup\n"
      "    structure of the assembly calculation to the few-group structure "
      "for\n"
      "    the nodal diffusion calculation. If None, no few-group diffusion\n"
      "    constants will be generated.\n"
      "pins : list of FuelPin or GuideTube\n"
      "    Description of the pins which make up the assembly geometry.\n"
      "dancoff_track_spacing : float\n"
      "    Spacing between tracks when calculating Dancoff factors.\n"
      "    Default is 0.05.\n"
      "dancoff_num_azimuthal_angles : int\n"
      "    Number of azimuthal angles when calculating Dancoff factors.\n"
      "    Default is 64.\n"
      "dancoff_polar_quadrature: PolarQuadrature\n"
      "    The polar quadrature used when calculating Dancoff factors.\n"
      "    Default is YamamotoTabuchi6.\n"
      "dancoff_isolation_factor : float\n"
      "    The factor used to multiply the pitch when calculating the flux in "
      "an\n"
      "    isolated pin. Default is 20.\n"
      "moderator_xs : CrossSection\n"
      "    The micro-group cross sections for the moderator in the problem.\n"
      "average_fuel_pin : CrossSection\n"
      "    Homogenized micro-group cross sections for the average fuel pin "
      "cell.\n"
      "track_spacing : float\n"
      "    Spacing between tracks in the assembly calculation. Default is "
      "0.02.\n"
      "num_azimuthal_angles : int\n"
      "    Number of azimuthal angles in the assembly calculation. Default is "
      "32.\n"
      "polar_quadrature: PolarQuadrature\n"
      "    The polar quadrature used in the assembly calculation.\n"
      "    Default is YamamotoTabuchi6.\n"
      "keff_tolerance : float\n"
      "    Convergence criteria for keff. Default is 1.E-5.\n"
      "flux_tolerance : float\n"
      "    Convergence criteria for the flux. Default is 1.E-5.\n"
      "plot_assembly : bool\n"
      "    Indicates wether the GUI plotter for the assembly geometry will be\n"
      "    activated before performing the calcualtion.\n"
      "fuel_dancoff_corrections : list of float\n"
      "    List of the dancoff corrections for the fuel in each pin.\n"
      "clad_dancoff_corrections : list of float\n"
      "    List of the dancoff corrections for the cladding in each pin.\n"
      "moc_geom : Cartesian2D\n"
      "    Geometry used in the assembly calculation. Is None until solve has\n"
      "    been called.\n"
      "moc : MOCDriver\n"
      "    The method of characteristics solver for the assembly calculation.\n"
      "    Is None until solve has been called.\n"
      "assmebly_diffusion_xs : DiffusionCrossSection\n"
      "    The few-group diffsuion group constants for the assembly.\n"
      "reflector_diffusion_xs : DiffusionCrossSection\n"
      "    The few-group diffsuion group constants for the reflector.\n"
      "adf : ndarray\n"
      "    The assembly discontinuity factors.\n"
      "cdf : ndarray\n"
      "    The corner discontinuity factors.\n")

      .def(py::init<double /*pitch*/, std::shared_ptr<Material> /*moderator*/,
                    std::pair<std::size_t, std::size_t> /*shape*/,
                    double /*gap_width*/, double /*baffle_width*/,
                    std::shared_ptr<Material> /*baffle*/,
                    std::shared_ptr<NDLibrary> /*ndl*/>(),
           py::arg("pitch"), py::arg("moderator"), py::arg("shape"),
           py::arg("gap_width"), py::arg("baffle_width"), py::arg("baffle"),
           py::arg("ndl"))

      .def("solve", &PWRReflector::solve,
           "Solve the assembly-reflector, generating few-group diffusion "
           "constants.")

      .def("save_diffusion_data", &PWRReflector::save_diffusion_data,
           "Saves the diffusion data to a numpy zip file.\n\n"
           "Parameters\n"
           "----------\n"
           "fname : str\n"
           "        Name of file in which to save data.",
           py::arg("fname"))

      .def_property("condensation_scheme", &PWRReflector::condensation_scheme,
                    &PWRReflector::set_condensation_scheme)

      .def_property("few_group_condensation_scheme",
                    &PWRReflector::few_group_condensation_scheme,
                    &PWRReflector::set_few_group_condensation_scheme)

      .def_property("pins", &PWRReflector::pins, &PWRReflector::set_pins)

      .def_property("moderator", &PWRReflector::moderator,
                    &PWRReflector::set_moderator)

      .def_property("num_azimuthal_angles", &PWRReflector::num_azimuthal_angles,
                    &PWRReflector::set_num_azimuthal_angles)

      .def_property("track_spacing", &PWRReflector::track_spacing,
                    &PWRReflector::set_track_spacing)

      .def_property("polar_quadrature", &PWRReflector::polar_quadrature,
                    &PWRReflector::set_polar_quadrature)

      .def_property("dancoff_num_azimuthal_angles",
                    &PWRReflector::dancoff_num_azimuthal_angles,
                    &PWRReflector::set_dancoff_num_azimuthal_angles)

      .def_property("dancoff_track_spacing",
                    &PWRReflector::dancoff_track_spacing,
                    &PWRReflector::set_dancoff_track_spacing)

      .def_property("dancoff_polar_quadrature",
                    &PWRReflector::dancoff_polar_quadrature,
                    &PWRReflector::set_dancoff_polar_quadrature)

      .def_property("keff_tolerance", &PWRReflector::keff_tolerance,
                    &PWRReflector::set_keff_tolerance)

      .def_property("flux_tolerance", &PWRReflector::flux_tolerance,
                    &PWRReflector::set_flux_tolerance)

      .def_property("plot_assembly", &PWRReflector::plot_assembly,
                    &PWRReflector::set_plot_assembly)

      .def_property_readonly("fuel_dancoff_corrections",
                             &PWRReflector::fuel_dancoff_corrections)

      .def_property_readonly("clad_dancoff_corrections",
                             &PWRReflector::clad_dancoff_corrections)

      .def_property_readonly("moderator_xs", &PWRReflector::moderator_xs)

      .def_property_readonly("average_fuel_pin",
                             &PWRReflector::average_fuel_pin)

      .def_property_readonly("adf", &PWRReflector::adf)

      .def_property_readonly("cdf", &PWRReflector::cdf)

      .def_property_readonly("assembly_diffusion_xs",
                             &PWRReflector::assembly_diffusion_xs)

      .def_property_readonly("reflector_diffusion_xs",
                             &PWRReflector::reflector_diffusion_xs)

      .def_property_readonly("moc", &PWRReflector::moc)

      .def_property_readonly("moc_geom", &PWRReflector::moc_geom)

      .def_property_readonly("pitch", &PWRReflector::pitch)

      .def_property_readonly("shape", &PWRReflector::pitch);
}
