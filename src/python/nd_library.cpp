#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/nd_library.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_NuclideHandle(py::module& m) {
  py::class_<NuclideHandle>(m, "NuclideHandle",
      "A NuclideHandle contains information about a specific nuclide, obtained "
      "from the nuclear data library (see :py:class:`NDLibrary` ).")

      .def_readonly("name", &NuclideHandle::name, "Indentifier of the nuclide.")

      .def_readonly("label", &NuclideHandle::label,
                    "Optional label provided at library creation.")

      .def_readonly(
          "temperatures", &NuclideHandle::temperatures,
          "List of temperatures at which cross sections are tabulated.")

      .def_readonly("dilutions", &NuclideHandle::dilutions,
                    "List of dilutions at which cross sections are tabulated.")

      .def_readonly("awr", &NuclideHandle::awr,
                    "Atomic weight ratio of the nuclide.")

      .def_readonly("potential_xs", &NuclideHandle::potential_xs,
                    "Potential scattering cross section of the nuclide.")

      .def_readonly("ZA", &NuclideHandle::ZA,
                    "The ZA number of the nuclide, constructed as Z*1000 + A.")

      .def_readonly("fissile", &NuclideHandle::fissile,
                    "True if the nuclide is fissile, False otherwise.")

      .def_readonly("resonant", &NuclideHandle::resonant,
                    "True if the nuclide is resonant, False otherwise.");
}

void init_NDLibrary(py::module& m) {
  py::class_<NDLibrary, std::shared_ptr<NDLibrary>>(m, "NDLibrary")
      .def(py::init<const std::string&>(),
           "Creates a new NDLibrary object from an HDF5 file.\n\n"
           "Parameters\n"
           "----------\n"
           "fname : str\n"
           "        Name of the hdf5 file with the library.",
           py::arg("fname"))

      .def("get_nuclide",
           py::overload_cast<const std::string&>(&NDLibrary::get_nuclide,
                                                 py::const_),
           "Returns the :py:class:`NuclideHandle` of the the desired nuclide.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the desired nuclide.",
           py::arg("name"))

      .def("interp_xs", &NDLibrary::interp_xs,
           "Interpolates the cross section of the prescribed nuclide to the "
           "desired temperature and dilution.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the desired nuclide.\n"
           "temp : str\n"
           "       Desired temperature in kelvin.\n"
           "dil : str\n"
           "      Desired dilution in barns.",
           py::arg("name"), py::arg("temp"), py::arg("dil"))

      .def("two_term_xs", &NDLibrary::two_term_xs,
           "Uses the two-term rational approximation for self shielding of "
           "cross sections, where the fuel escape probability is approximated "
           "as \n\n"
           ".. math:: P^{F\\to M}(E) = \\beta_1\\frac{\\alpha_1 \\Sigma_e}{\\Sigma_t(E) + \\alpha_1\\Sigma_e} + "
           "\\beta_2\\frac{\\alpha_2 \\Sigma_e}{\\Sigma_t(E) + \\alpha_2\\Sigma_e}.\n\n"
           "If shielding isotope :math:`r` with potential cross section "
           ":math:`\\sigma_p^r` and number density :math:`N_r`, then the "
           "background cross sections are computed as\n\n"
           ".. math:: \\sigma_1 = \\frac{\\Sigma_p - N_r\\sigma_p^r + \\alpha_1\\Sigma_e}{N_r} \\\\\\\\"
           "          \\sigma_2 = \\frac{\\Sigma_p - N_r\\sigma_p^r + \\alpha_2\\Sigma_e}{N_r}\n\n"
           "where :math:`\\Sigma_p` is the macroscopic potential scattering "
           "cross section of the material and :math:`\\Sigma_e` is the escape "
           "cross section.\n\n"
           "Parameters\n"
           "----------\n"
           "name : str\n"
           "       Name of the nuclide to be treated.\n"
           "temp : float\n"
           "       Temperature of the material (in kelvin).\n"
           "a1 : float\n"
           "     :math:`\\alpha_1`.\n"
           "a2 : float\n"
           "     :math:`\\alpha_2`.\n"
           "b1 : float\n"
           "     :math:`\\beta_1`.\n"
           "b2 : float\n"
           "     :math:`\\beta_2`.\n"
           "xs1 : float\n"
           "      Background cross section for first term (:math:`\\sigma_1`).\n"
           "xs2 : float\n"
           "      Background cross section for second term (:math:`\\sigma_2`).",
           py::arg("name"), py::arg("temp"), py::arg("a1"), py::arg("a2"),
           py::arg("b1"), py::arg("b2"), py::arg("xs1"), py::arg("xs2"))

      .def_property_readonly("library", &NDLibrary::library,
                             "Name of the nuclear data library (if provided).")

      .def_property_readonly("ngroups", &NDLibrary::ngroups,
                             "Number of energy groups in the library.")

      .def_property_readonly("group_bounds", &NDLibrary::group_bounds,
                             "The boundaries of the energy groups for the "
                             "group structure (in decreasing order).")

      .def_property_readonly("group_structure", &NDLibrary::group_structure,
                             "The name of the group structure (if provided).");
}