#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <utils/nuclide_names.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_NuclideNameFuncs(py::module& m) {
  m.def("nuclide_name_to_simple_name", &nuclide_name_to_simple_name,
        "Removes any scattering law names appended to nuclide. For example, \n"
        "H1_H2O becomes H1.\n\n"
        "Parameters\n"
        "----------\n"
        "name : string\n"
        "    Name of the nuclide.\n\n"
        "Returns\n"
        "-------\n"
        "string\n"
        "    Simplified name of the nuclide.\n",
        py::arg("name"));

  m.def("nuclide_name_to_element_symbol", &nuclide_name_to_element_symbol,
        "Returns only the element symbol of a nuclide name. For example, if "
        "provided with \"Ag109\", only \"Ag\" will be returned.\n\n"
        "Parameters\n"
        "----------\n"
        "name : string\n"
        "    Name of the nuclide.\n\n"
        "Returns\n"
        "-------\n"
        "string\n"
        "    Element symbol of the nuclide.\n",
        py::arg("name"));

  m.def("nuclide_name_to_za", &nuclide_name_to_za,
        "Returns the ZAID of a nuclide. For example, if provided with "
        "\"U235\", 92235 will be returned. This also provides a unique ZAID "
        "for isomers such as Am242m1.\n\n"
        "Parameters\n"
        "----------\n"
        "name : string\n"
        "    Name of the nuclide.\n\n"
        "Returns\n"
        "-------\n"
        "int\n"
        "    ZAID of the nuclide.\n",
        py::arg("name"));
}