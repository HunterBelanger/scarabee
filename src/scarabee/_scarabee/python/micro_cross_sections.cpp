#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <data/micro_cross_sections.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_MicroCrossSectionStructs(py::module& m) {
  py::class_<MicroDepletionXS>(
      m, "MicroDepletionXS",
      "Holds the microscopic cross sections needed for depletion calculations.")
      .def_readonly("n_fission", &MicroDepletionXS::n_fission,
                    "Fission cross section.")
      .def_readonly("n_gamma", &MicroDepletionXS::n_gamma,
                    "Capture cross section.")
      .def_readonly("n_2n", &MicroDepletionXS::n_2n, "(n,2n) cross section.")
      .def_readonly("n_3n", &MicroDepletionXS::n_3n, "(n,3n) cross section.")
      .def_readonly("n_a", &MicroDepletionXS::n_a, "(n,a) cross section.")
      .def_readonly("n_p", &MicroDepletionXS::n_p, "(n,p) cross section.");

  py::class_<MicroNuclideXS>(
      m, "MicroNuclideXS",
      "Holds the microscopic cross sections for a nuclide.")
      .def_readonly("Et", &MicroNuclideXS::Et, "Total cross section.")
      .def_readonly("Dtr", &MicroNuclideXS::Dtr, "Transport correction.")
      .def_readonly("Es", &MicroNuclideXS::Es, "Scattering cross section.")
      .def_readonly("Ea", &MicroNuclideXS::Ea, "Absorption cross section.")
      .def_readonly("Ef", &MicroNuclideXS::Ef, "Fission cross section.")
      .def_readonly("nu", &MicroNuclideXS::nu, "Fission yield.")
      .def_readonly("chi", &MicroNuclideXS::chi, "Fission spectrum.");

  py::class_<ResonantOneGroupXS>(m, "ResonantOneGroupXS",
                                 "Single group resonant cross sections.")
      .def_readonly("Dtr", &ResonantOneGroupXS::Dtr, "Transport correction.")
      .def_readonly("Ea", &ResonantOneGroupXS::Ea, "Absorption cross section.")
      .def_readonly("Ef", &ResonantOneGroupXS::Ef, "Fission cross section.")
      .def_readonly("Es", &ResonantOneGroupXS::Es, "Scattering cross section.")
      .def_readonly("gout_min", &ResonantOneGroupXS::gout_min,
                    "First outgoing energy group.")
      .def_readonly("n_gamma", &ResonantOneGroupXS::n_gamma,
                    "Capture cross section.");
}
