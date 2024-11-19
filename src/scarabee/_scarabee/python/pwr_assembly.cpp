#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <assemblies/pwr_assembly.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_PWRAssembly(py::module& m) {
  py::class_<PWRAssembly, std::shared_ptr<PWRAssembly>>(
      m, "PWRAssembly",
      "A PWRAssembly instance is responsible for performing all the lattice\n"
      "calculations necessary to produce few-group cross sections for a "
      "single\n"
      "PWR assembly.\n\n"
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
      "boundary_conditions : BoundaryCondition\n"
      "    Boundary condition to be applied to all sides of the assembly.\n"
      "    Default value is Periodic.\n"
      "plot_assembly : bool\n"
      "    Indicates wether the GUI plotter for the assembly geometry will be\n"
      "    activated before performing the calcualtion.\n"
      "fuel_dancoff_corrections : list of float\n"
      "    List of the dancoff corrections for the fuel in each pin.\n"
      "clad_dancoff_corrections : list of float\n"
      "    List of the dancoff corrections for the cladding in each pin.\n"
      "criticality_spectrum_method : str or None\n"
      "    The type of leakage approximation used to modify the assembly flux\n"
      "    spectrum. Acceptable values are \"B1\", \"P1\", or None.\n"
      "moc_geom : Cartesian2D\n"
      "    Geometry used in the assembly calculation. Is None until solve has\n"
      "    been called.\n"
      "moc : MOCDriver\n"
      "    The method of characteristics solver for the assembly calculation.\n"
      "    Is None until solve has been called.\n"
      "diffusion_xs : DiffusionCrossSection\n"
      "    The few-group diffsuion group constants, if few-group constants "
      "were\n"
      "    generated.\n"
      "form_factors : ndarray\n"
      "    The pin-wise form factors for power distribution reconstruction.\n"
      "    Is None until solve has been called.\n"
      "adf : ndarray\n"
      "    The assembly discontinuity factors.\n"
      "cdf : ndarray\n"
      "    The corner discontinuity factors.\n")

      .def(py::init<double /*pitch*/, std::shared_ptr<Material> /*moderator*/,
                    std::pair<std::size_t, std::size_t> /*shape*/,
                    std::shared_ptr<NDLibrary> /*ndl*/>(),
           py::arg("pitch"), py::arg("moderator"), py::arg("shape"),
           py::arg("ndl"))

      .def("solve", &PWRAssembly::solve,
           "Solve the assembly, generating few-group diffusion constants.")

      .def("save_diffusion_data", &PWRAssembly::save_diffusion_data,
           "Saves the diffusion data to a numpy zip file.\n\n"
           "Parameters\n"
           "----------\n"
           "fname : str\n"
           "        Name of file in which to save data.",
           py::arg("fname"))

      .def_property("criticality_spectrum_method",
                    &PWRAssembly::criticality_spectrum_method,
                    &PWRAssembly::set_criticality_spectrum_method)

      .def_property("condensation_scheme", &PWRAssembly::condensation_scheme,
                    &PWRAssembly::set_condensation_scheme)

      .def_property("few_group_condensation_scheme",
                    &PWRAssembly::few_group_condensation_scheme,
                    &PWRAssembly::set_few_group_condensation_scheme)

      .def_property("pins", &PWRAssembly::pins, &PWRAssembly::set_pins)

      .def_property("moderator", &PWRAssembly::moderator,
                    &PWRAssembly::set_moderator)

      .def_property("num_azimuthal_angles", &PWRAssembly::num_azimuthal_angles,
                    &PWRAssembly::set_num_azimuthal_angles)

      .def_property("track_spacing", &PWRAssembly::track_spacing,
                    &PWRAssembly::set_track_spacing)

      .def_property("polar_quadrature", &PWRAssembly::polar_quadrature,
                    &PWRAssembly::set_polar_quadrature)

      .def_property("dancoff_num_azimuthal_angles",
                    &PWRAssembly::dancoff_num_azimuthal_angles,
                    &PWRAssembly::set_dancoff_num_azimuthal_angles)

      .def_property("dancoff_track_spacing",
                    &PWRAssembly::dancoff_track_spacing,
                    &PWRAssembly::set_dancoff_track_spacing)

      .def_property("dancoff_polar_quadrature",
                    &PWRAssembly::dancoff_polar_quadrature,
                    &PWRAssembly::set_dancoff_polar_quadrature)

      .def_property("keff_tolerance", &PWRAssembly::keff_tolerance,
                    &PWRAssembly::set_keff_tolerance)

      .def_property("flux_tolerance", &PWRAssembly::flux_tolerance,
                    &PWRAssembly::set_flux_tolerance)

      .def_property("boundary_conditions", &PWRAssembly::boundary_conditions,
                    &PWRAssembly::set_boundary_conditions)

      .def_property("plot_assembly", &PWRAssembly::plot_assembly,
                    &PWRAssembly::set_plot_assembly)

      .def_property_readonly("fuel_dancoff_corrections",
                             &PWRAssembly::fuel_dancoff_corrections)

      .def_property_readonly("clad_dancoff_corrections",
                             &PWRAssembly::clad_dancoff_corrections)

      .def_property_readonly("moderator_xs", &PWRAssembly::moderator_xs)

      .def_property_readonly("average_fuel_pin", &PWRAssembly::average_fuel_pin)

      .def_property_readonly("form_factors", &PWRAssembly::form_factors)

      .def_property_readonly("adf", &PWRAssembly::adf)

      .def_property_readonly("cdf", &PWRAssembly::cdf)

      .def_property_readonly("diffusion_xs", &PWRAssembly::diffusion_xs)

      .def_property_readonly("moc", &PWRAssembly::moc)

      .def_property_readonly("moc_geom", &PWRAssembly::moc_geom)

      .def_property_readonly("pitch", &PWRAssembly::pitch)

      .def_property_readonly("shape", &PWRAssembly::pitch);
}
