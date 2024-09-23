#include <pybind11/pybind11.h>

#include <moc/track.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_Track(py::module& m) {
  py::class_<Track>(m, "Track")
      .def("wgt", &Track::wgt)
      .def("width", &Track::width)
      .def("phi", &Track::phi)
      .def("entry_pos", &Track::entry_pos)
      .def("exit_pos", &Track::exit_pos)
      .def("entry_flux", py::overload_cast<>(&Track::entry_flux, py::const_))
      .def("exit_flux", py::overload_cast<>(&Track::exit_flux, py::const_));
}
