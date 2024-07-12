#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>
#include <ImApp/imapp.hpp>

#include <moc/moc_driver.hpp>
#include <moc/moc_plotter.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <tuple>

namespace py = pybind11;

using namespace scarabee;

void init_MOCDriver(py::module& m) {
  py::class_<MOCDriver>(m, "MOCDriver")
      .def(py::init<
               std::shared_ptr<Cartesian2D> /*geometry*/,
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
           "    ymaxbc    BoundaryCondition for maximum y value",
           py::arg("geometry"),
           py::arg("xminbc") = BoundaryCondition::Reflective,
           py::arg("xmaxbc") = BoundaryCondition::Reflective,
           py::arg("yminbc") = BoundaryCondition::Reflective,
           py::arg("ymaxbc") = BoundaryCondition::Reflective)

      .def("generate_tracks",
           py::overload_cast<std::uint32_t, double, PolarQuadrature>(
               &MOCDriver::generate_tracks),
           "Traces tracks for the calculation across the geometry.\n\n"
           "Arguments:\n"
           "    nangles       number of azimuthal angles (even)\n"
           "    d             max spacing between tracks of a given angle\n"
           "    polar_quad    polar quadrature for generating segment lengths",
           py::arg("nangles"), py::arg("d"), py::arg("polar_quad"))

      .def(
          "generate_tracks",
          [](MOCDriver& md, std::uint32_t na, double d,
             PolarQuadratureType pq) { return md.generate_tracks(na, d, pq); },
          "Traces tracks for the calculation across the geometry.\n\n"
          "Arguments:\n"
          "    nangles       number of azimuthal angles (even)\n"
          "    d             max spacing between tracks of a given angle\n"
          "    polar_quad    polar quadrature for generating segment lengths",
          py::arg("nangles"), py::arg("d"), py::arg("polar_quad"))

      .def("drawn", &MOCDriver::drawn,
           "return True if geometry has been traced")

      .def_property(
          "keff_tolerance", &MOCDriver::keff_tolerance,
          &MOCDriver::set_keff_tolerance,
          "maximum relative absolute difference in keff for convergence")

      .def_property(
          "flux_tolerance", &MOCDriver::flux_tolerance,
          &MOCDriver::set_flux_tolerance,
          "maximum relative absolute difference in flux for convergence")

      .def("flux",
           py::overload_cast<const Vector&, const Direction&, std::size_t>(
               &MOCDriver::flux, py::const_),
           "Returns the scalar flux in group g at position r.\n\n"
           "Arguments:\n"
           "    r  position\n"
           "    u  direction\n"
           "    g  group index",
           py::arg("r"), py::arg("u"), py::arg("g"))

      .def("flux",
           py::overload_cast<std::size_t, std::size_t>(&MOCDriver::flux,
                                                       py::const_),
           "Returns the scalar flux in group g in FSR i.\n\n"
           "Arguments:\n"
           "    i  FSR index\n"
           "    g  group index",
           py::arg("i"), py::arg("g"))

      .def("volume",
           py::overload_cast<const Vector&, const Direction&>(
               &MOCDriver::volume, py::const_),
           "Returns the volume of the flat source region at position r.\n\n"
           "Arguments:\n"
           "    r  position\n"
           "    u  direction",
           py::arg("r"), py::arg("u"))

      .def("volume",
           py::overload_cast<std::size_t>(&MOCDriver::volume, py::const_),
           "Returns the volume of Flat Source Region i.\n\n"
           "Arguments:\n"
           "    i  FSR index",
           py::arg("i"))

      .def("xs",
           py::overload_cast<const Vector&, const Direction&>(&MOCDriver::xs,
                                                              py::const_),
           "Returns the CrossSection at position r.\n\n"
           "Arguments:\n"
           "    r  position\n"
           "    u  direction",
           py::arg("r"), py::arg("u"))

      .def("xs", py::overload_cast<std::size_t>(&MOCDriver::xs, py::const_),
           "Returns the CrossSection in FSR i.\n\n"
           "Arguments:\n"
           "    i  FSR index",
           py::arg("i"))

      .def_property_readonly("keff", &MOCDriver::keff,
                             "value of keff  estimated by solver (1 by default "
                             "if no solution has been obtained)")

      .def_property_readonly("ngroups", &MOCDriver::ngroups,
                             "number of energy groups")

      .def_property_readonly("polar_quadrature", &MOCDriver::polar_quadrature,
                             "quadrature used for polar angle integration")

      .def("solve", &MOCDriver::solve, "begins iterations to solve problem")

      .def_property(
          "sim_mode",
          [](const MOCDriver& md) -> SimulationMode { return md.sim_mode(); },
          [](MOCDriver& md, SimulationMode& m) { md.sim_mode() = m; },
          "simulation mode (fixed-source, keff)")

      .def_property(
          "x_min_bc",
          [](const MOCDriver& md) -> BoundaryCondition {
            return md.x_min_bc();
          },
          [](MOCDriver& md, BoundaryCondition& bc) { md.x_min_bc() = bc; },
          "boundadary condition at x_min")

      .def_property(
          "x_max_bc",
          [](const MOCDriver& md) -> BoundaryCondition {
            return md.x_max_bc();
          },
          [](MOCDriver& md, BoundaryCondition& bc) { md.x_max_bc() = bc; },
          "boundadary condition at x_max")

      .def_property(
          "y_min_bc",
          [](const MOCDriver& md) -> BoundaryCondition {
            return md.y_min_bc();
          },
          [](MOCDriver& md, BoundaryCondition& bc) { md.y_min_bc() = bc; },
          "boundadary condition at y_min")

      .def_property(
          "y_max_bc",
          [](const MOCDriver& md) -> BoundaryCondition {
            return md.y_max_bc();
          },
          [](MOCDriver& md, BoundaryCondition& bc) { md.y_max_bc() = bc; },
          "boundadary condition at y_max")

      .def_property_readonly("geometry", &MOCDriver::geometry,
                             "the Cartesian2D geometry for the problem")

      .def_property_readonly("size", &MOCDriver::size,
                             "number of flat source regions")

      .def_property_readonly("nfsr", &MOCDriver::nfsr,
                             "number of flat source regions")

      .def_property_readonly("nregions", &MOCDriver::nregions,
                             "number of flat source regions")

      .def_property_readonly("x_min", &MOCDriver::x_min,
                             "minimum value of x in problem domain")

      .def_property_readonly("x_max", &MOCDriver::x_max,
                             "maximum value of x in problem domain")

      .def_property_readonly("y_min", &MOCDriver::y_min,
                             "minimum value of y in problem domain")

      .def_property_readonly("y_max", &MOCDriver::y_max,
                             "maximum value of y in problem domain")

      .def_property_readonly("solved", &MOCDriver::solved,
                             "returns true if solve_keff has been run "
                             "sucessfully (reset to false on generate_tracks)")

      .def("set_extern_src",
           py::overload_cast<const Vector&, const Direction&, std::size_t,
                             double>(&MOCDriver::set_extern_src),
           "Sets the external source in the flat source region at r.\n\n"
           "Arguments:\n"
           "    r    position\n"
           "    u    direction\n"
           "    g    group index\n"
           "    src  value of source",
           py::arg("r"), py::arg("u"), py::arg("g"), py::arg("src"))

      .def("set_extern_src",
           py::overload_cast<std::size_t, std::size_t, double>(
               &MOCDriver::set_extern_src),
           "Sets the external source in flat source region i.\n\n"
           "Arguments:\n"
           "    i    flat source region index\n"
           "    g    group index\n"
           "    src  value of source",
           py::arg("i"), py::arg("g"), py::arg("src"))

      .def("extern_src",
           py::overload_cast<const Vector&, const Direction&, std::size_t>(
               &MOCDriver::extern_src, py::const_),
           "Returns the external source in the flat source region at r.\n\n"
           "Arguments:\n"
           "    r    position\n"
           "    u    direction\n"
           "    g    group index",
           py::arg("r"), py::arg("u"), py::arg("g"))

      .def("extern_src",
           py::overload_cast<std::size_t, std::size_t>(&MOCDriver::extern_src,
                                                       py::const_),
           "Returns the external source in flat source region i.\n\n"
           "Arguments:\n"
           "    i    flat source region index\n"
           "    g    group index",
           py::arg("i"), py::arg("g"))

      .def("homogenize",
           py::overload_cast<>(&MOCDriver::homogenize, py::const_),
           "Computes a homogenized set of cross sections for the problem.\n\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "            Homogenized cross section.")

      .def("homogenize",
           py::overload_cast<const std::vector<std::size_t>&>(
               &MOCDriver::homogenize, py::const_),
           "Computes a homogenized set of cross sections for all provided "
           "regions.\n\n"
           "Parameters\n"
           "----------\n"
           "regions : list of int\n"
           "        List of regions for homogenization.\n"
           "Returns\n"
           "-------\n"
           "CrossSection\n"
           "            Homogenized cross section.",
           py::arg("regions"))

      .def(
          "homogenize_flux_spectrum",
          py::overload_cast<>(&MOCDriver::homogenize_flux_spectrum, py::const_),
          "Computes a homogenized flux spectrum which can be used for energy "
          "condensation.\n\n"
          "Returns\n"
          "-------\n"
          "ndarray of floats\n"
          "                 Homogenized flux spectrum.")

      .def("homogenize_flux_spectrum",
           py::overload_cast<const std::vector<std::size_t>&>(
               &MOCDriver::homogenize_flux_spectrum, py::const_),
           "Computes a homogenized flux spectrum which can be used for energy "
           "condensation for all provided regions.\n\n"
           "Parameters\n"
           "----------\n"
           "regions : list of int\n"
           "        List of regions for homogenization.\n"
           "Returns\n"
           "-------\n"
           "ndarray of floats\n"
           "                 Homogenized flux spectrum.",
           py::arg("regions"))

      .def("apply_criticality_spectrum", &MOCDriver::apply_criticality_spectrum,
           "Modifies the flux spectrum of the solved problem by multiplying "
           "the value of the flux by the ratio of the provided criticality "
           "spectrum to homogenized flux spectrum.\n\n"
           "Parameters\n"
           "----------\n"
           "flux : ndarray of floats\n"
           "       Criticality spectrum from a P1 or B1 calculation.\n",
           py::arg("flux"))

      .def("plot",
           [](const MOCDriver& md) {
             ImApp::App guiplotter(1920, 1080, "Scarabee MOC Plotter");
             guiplotter.enable_docking();
             guiplotter.push_layer(std::make_unique<MOCPlotter>(&md));
             guiplotter.run();
           })

      .def(
          "rasterize_flux",
          [](const MOCDriver& md, std::size_t nx, std::size_t ny) {
            if (nx == 0) {
              auto mssg = "Rasterizing flux requires nx > 0.";
              spdlog::error(mssg);
              throw ScarabeeException(mssg);
            } else if (ny == 0) {
              auto mssg = "Rasterizing flux requires ny > 0.";
              spdlog::error(mssg);
              throw ScarabeeException(mssg);
            }

            if (md.solved() == false) {
              auto mssg = "Cannot rasterize flux. System has not been solved.";
              spdlog::error(mssg);
              throw ScarabeeException(mssg);
            }

            const double dx =
                (md.x_max() - md.x_min()) / static_cast<double>(nx);
            const double dy =
                (md.y_max() - md.y_min()) / static_cast<double>(ny);
            const Direction u(1., 1.);

            // Initialize arrays to return
            xt::pytensor<double, 3> flux =
                xt::zeros<double>({md.ngroups(), ny, nx});
            xt::pytensor<double, 1> x_bnds = xt::zeros<double>({nx + 1});
            xt::pytensor<double, 1> y_bnds = xt::zeros<double>({ny + 1});

            // Fill bounds
            x_bnds(0) = md.x_min();
            for (std::size_t i = 0; i < nx; i++) {
              x_bnds(i + 1) = x_bnds(i) + dx;
            }

            y_bnds(0) = md.y_min();
            for (std::size_t j = 0; j < ny; j++) {
              y_bnds(j + 1) = y_bnds(j) + dy;
            }

            // Rasterize flux
            for (std::size_t g = 0; g < md.ngroups(); g++) {
              double y = md.y_min() + 0.5 * dy;
              for (std::size_t j = 0; j < ny; j++) {
                double x = md.x_min() + 0.5 * dx;
                for (std::size_t i = 0; i < nx; i++) {
                  Vector r(x, y);
                  flux(g, j, i) = md.flux(r, u, g);
                  x += dx;
                }
                y += dy;
              }
            }

            // Return values
            return std::tuple<xt::pytensor<double, 3>, xt::pytensor<double, 1>,
                              xt::pytensor<double, 1>>{flux, x_bnds, y_bnds};
          },
          "Rasterizes the flux in all energy groups for easy plotting.\n\n"
          "Arguments:\n"
          "    nx  number of mesh bins along x\n"
          "    ny  number of mesh bins along y");
}
