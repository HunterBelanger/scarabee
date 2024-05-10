#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/nd_library.hpp>

#include <memory>

namespace py = pybind11;

using namespace scarabee;

void init_NuclideHandle(py::module& m) {
  py::class_<NuclideHandle>(m, "NuclideHandle")
  .def_readonly("name", &NuclideHandle::name)
  .def_readonly("label", &NuclideHandle::label)
  .def_readonly("temperatures", &NuclideHandle::temperatures)
  .def_readonly("dilutions", &NuclideHandle::dilutions)
  .def_readonly("awr", &NuclideHandle::awr)
  .def_readonly("potential_xs", &NuclideHandle::potential_xs)
  .def_readonly("ZA", &NuclideHandle::ZA)
  .def_readonly("fissile", &NuclideHandle::fissile)
  .def_readonly("resonant", &NuclideHandle::resonant);
}

void init_NDLibrary(py::module& m) {
  py::class_<NDLibrary, std::shared_ptr<NDLibrary>>(m, "NDLibrary")
  .def(py::init<const std::string&>())
  .def_property_readonly("library", &NDLibrary::library)
  .def_property_readonly("ngroups", &NDLibrary::ngroups)
  .def_property_readonly("group_structure", &NDLibrary::group_structure)
  .def("get_nuclide", py::overload_cast<const std::string&>(&NDLibrary::get_nuclide, py::const_))
  .def("interp_nuclide_xs", &NDLibrary::interp_nuclide_xs)
  .def("carlvik_two_term", &NDLibrary::carlvik_two_term)
  .def("potential_xs", &NDLibrary::potential_xs);
}